// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EIGEN_STEPPER_HPP
#define ACTS_EIGEN_STEPPER_HPP 1

#include <cmath>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/concept/AnyFieldLookup.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

ActsMatrixD<3, 3>
cross(const ActsMatrixD<3, 3>& m, const Vector3D& v)
{
  ActsMatrixD<3, 3> r;
  r.col(0) = m.col(0).cross(v);
  r.col(1) = m.col(1).cross(v);
  r.col(2) = m.col(2).cross(v);

  return r;
}

/// Runge-Kutta-Nystroem stepper for the following ODE:
///
/// r = (x,y,z)    ... global position
/// T = (Ax,Ay,Az) ... momentum direction (normalized)
///
/// dr/ds = T
/// dT/ds = q/p * (T x B)
///
/// with s being the arc length of the track, q the charge of the particle,
/// p its momentum and B the magnetic field
///
template <typename BField>
class EigenStepper
{

private:
  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    typedef BoundParameters type;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    typedef CurvilinearParameters type;
  };

  /// Rotation matrix going from global coordinates to local surface coordinates
  static ActsMatrixD<3, 3>
  dLocaldGlobal(const Surface& p, const Vector3D& gpos)
  {
    // A surface's associated transform is a 4x4 affine transformation matrix,
    // whose top-left corner is the 3x3 linear local-to-global rotation matrix,
    // and whose right column is a translation vector, with a trailing 1.
    //
    // So to get the global-to-local rotation matrix, we only need to take the
    // transpose of the top-left corner of the surface's associated transform.
    //
    return p.transform().matrix().topLeftCorner<3, 3>().transpose();
  }

public:
  /// Cache for track parameter propagation
  ///
  /// it is exposed to public for use of the expert-only
  /// propagate_with_cache method of the propagator
  struct Cache
  {
    /// Constructor from the initial track parameters
    /// @tparam [in] par The track parameters at start
    ///
    /// @note the covariance matrix is copied when needed
    template <typename T>
    explicit Cache(const T& par)
      : pos(par.position())
      , dir(par.momentum().normalized())
      , qop(par.charge() / par.momentum().norm())
      , cov_transport(false)
      , accumulated_path(0.)
      , step_size(std::numeric_limits<double>::max())
    {
      update_covariance(par);
    }

    /// The cache update for optimal performance
    /// @tparam [in] par The new track parameters at start
    ///
    /// @todo check to identify an reuse of start/cache
    template <typename T>
    void
    update(const T& par)
    {
      pos           = par.position();
      dir           = par.momentum().normalized();
      qop           = (par.charge() / par.momentum().norm());
      cov_transport = false;
      update_covariance(par);
    }

    /// The covariance update
    /// @tparam [in] par The (new) track parameters at start
    template <typename T>
    void
    update_covariance(const T& par)
    {
      if (par.covariance()) {
        cov_transport = true;
        cov           = ActsSymMatrixD<5>(*par.covariance());

        // The trigonometry required to convert the direction to spherical
        // coordinates and then compute the sines and cosines again can be
        // surprisingly expensive from a performance point of view.
        //
        // Here, we can avoid it because the direction is by definition a unit
        // vector, with the following coordinate conversions...
        //
        const double x = dir(0);  // == cos(phi) * sin(theta)
        const double y = dir(1);  // == sin(phi) * sin(theta)
        const double z = dir(2);  // == cos(theta)
        //
        // ...which we can invert to directly get the sines and cosines we want:
        //
        const double cos_theta     = z;
        const double sin_theta     = sqrt(x * x + y * y);
        const double inv_sin_theta = 1. / sin_theta;
        const double cos_phi       = x * inv_sin_theta;
        const double sin_phi       = y * inv_sin_theta;

        // @todo - check this might have to be the measurement frame ...
        const auto transform = par.referenceSurface().transform().matrix();
        jacobian(0, eLOC_0) = transform(0, eLOC_0);
        jacobian(0, eLOC_1) = transform(0, eLOC_1);
        jacobian(1, eLOC_0) = transform(1, eLOC_0);
        jacobian(1, eLOC_1) = transform(1, eLOC_1);
        jacobian(2, eLOC_0) = transform(2, eLOC_0);
        jacobian(2, eLOC_1) = transform(2, eLOC_1);
        jacobian(3, ePHI)   = -sin_theta * sin_phi;
        jacobian(3, eTHETA) = cos_theta * cos_phi;
        jacobian(4, ePHI)   = sin_theta * cos_phi;
        jacobian(4, eTHETA) = cos_theta * sin_phi;
        jacobian(5, eTHETA) = -sin_theta;
        jacobian(6, eQOP)   = 1;
      }
    }

    /// Global particle position accessor
    Vector3D
    position() const
    {
      return pos;
    }

    /// Momentum direction accessor
    Vector3D
    direction() const
    {
      return dir;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);

    /// Charge-momentum ratio, in natural units
    double qop = 1;

    /// Jacobian used to transport the covariance matrix
    ActsMatrixD<7, 5> jacobian = ActsMatrixD<7, 5>::Zero();

    /// The propagation derivative
    ActsVectorD<7> derivative = ActsVectorD<7>::Zero();

    /// Covariance matrix (and indicator)
    //// assocated with the initial error on track parameters
    bool              cov_transport = false;
    ActsSymMatrixD<5> cov           = ActsSymMatrixD<5>::Zero();

    /// Lazily initialized cache
    /// It caches the current magneticl field cell and stays interpolates within
    /// as long as this is valid. See step() code for details.
    bool                    field_cache_ready = false;
    concept::AnyFieldCell<> field_cache;

    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive sep size of the runge-kutta integration
    double step_size = std::numeric_limits<double>::max();
  };

  /// Always use the same propagation cache type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using cache_type = Cache;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  /// Return parameter types depend on the propagation mode:
  /// -  when propagating to a surface we usually return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor requires knowledge of the detector's magnetic field
  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  /// Convert the propagation cache to curvilinear parameters
  /// @param cache is the stepper cache
  /// @todo check: what if cache is already in curvilinear ? is this caught ?
  static CurvilinearParameters
  convert(const Cache& cache)
  {
    double                                   charge = cache.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov    = nullptr;

    // Perform error propagation if an initial covariance matrix was provided
    if (cache.cov_transport) {
      // Optimized trigonometry on the propagation direction, see documentation
      // of Cache::update_covariance() for a longer mathematical discussion.
      const double x = cache.dir(0);  // == cos(phi) * sin(theta)
      const double y = cache.dir(1);  // == sin(phi) * sin(theta)
      const double z = cache.dir(2);  // == cos(theta)
      //
      const double cos_theta     = z;
      const double sin_theta     = sqrt(x * x + y * y);
      const double inv_sin_theta = 1. / sin_theta;
      const double cos_phi       = x * inv_sin_theta;
      const double sin_phi       = y * inv_sin_theta;

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cos_theta) < 0.99) {
        // We normally operate in curvilinear coordinates defined as follows
        J(0, 0) = -sin_phi;
        J(0, 1) = cos_phi;
        J(1, 0) = -cos_phi * cos_theta;
        J(1, 1) = -sin_phi * cos_theta;
        J(1, 2) = sin_theta;
      } else {
        // Under grazing incidence to z, the above coordinate system definition
        // becomes numerically unstable, and we need to switch to another one
        const double c     = sqrt(y * y + z * z);
        const double inv_c = 1. / c;
        J(0, 1) = -z * inv_c;
        J(0, 2) = y * inv_c;
        J(1, 0) = c;
        J(1, 1) = -x * y * inv_c;
        J(1, 2) = -x * z * inv_c;
      }
      J(2, 3) = -sin_phi * inv_sin_theta;
      J(2, 4) = cos_phi * inv_sin_theta;
      J(3, 5) = -inv_sin_theta;
      J(4, 6) = 1;

      const ActsRowVectorD<5> scale_factors = cache.dir.transpose()
          * cache.jacobian.template topLeftCorner<3, 5>();

      const ActsMatrixD<5, 5> jac
          = J * (cache.jacobian - cache.derivative * scale_factors);

      cov = std::make_unique<const ActsSymMatrixD<5>>(jac * cache.cov
                                                      * jac.transpose());
    }

    return CurvilinearParameters(
        std::move(cov), cache.pos, cache.dir / std::abs(cache.qop), charge);
  }

  /// Convert the propagation cache to track parameters at a certain surface
  /// @param [in] cache Propagation cache used
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename S>
  static BoundParameters
  convert(Cache& cache, const S& surface)
  {
    double                                   charge = cache.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov    = nullptr;

    // Perform error propagation if an initial covariance matrix was provided
    if (cache.cov_transport) {
      // Optimized trigonometry on the propagation direction, see documentation
      // of Cache::update_covariance() for a longer mathematical discussion.
      const double x = cache.dir(0);  // == cos(phi) * sin(theta)
      const double y = cache.dir(1);  // == sin(phi) * sin(theta)
      //
      const double inv_sin_theta_2        = 1. / (x * x + y * y);
      const double cos_phi_over_sin_theta = x * inv_sin_theta_2;
      const double sin_phi_over_sin_theta = y * inv_sin_theta_2;
      const double inv_sin_theta          = sqrt(inv_sin_theta_2);

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      const auto dLdG = dLocaldGlobal(surface, cache.pos);
      J.topLeftCorner<2, 3>() = dLdG.template topLeftCorner<2, 3>();
      J(2, 3)                 = -sin_phi_over_sin_theta;
      J(2, 4)                 = cos_phi_over_sin_theta;
      J(3, 5)                 = -inv_sin_theta;
      J(4, 6)                 = 1;

      ActsRowVectorD<3> norm_vec = dLdG.template block<1, 3>(2, 0);
      norm_vec /= (norm_vec * cache.dir);

      const ActsRowVectorD<5> scale_factors
          = norm_vec * cache.jacobian.template topLeftCorner<3, 5>();

      const ActsMatrixD<5, 5> jac
          = J * (cache.jacobian - cache.derivative * scale_factors);

      cov = std::make_unique<const ActsSymMatrixD<5>>(jac * cache.cov
                                                      * jac.transpose());
    }

    return BoundParameters(std::move(cov),
                           cache.pos,
                           cache.dir / std::abs(cache.qop),
                           charge,
                           surface);
  }

  /// Estimate the (signed) distance to a certain surface
  static double
  distance(const Surface& s, const Vector3D& pos, const Vector3D& dir)
  {
    return s.intersectionEstimate(pos, dir).pathLength;
  }

  /// Get the field for the stepping
  /// It checks first if the access is still within the Cell,
  /// and updates the cell if necessary, then it takes the field
  /// from the cell
  /// @param [in,out] cache is the propagation cache associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(Cache& cache, const Vector3D& pos) const
  {
    if (!cache.field_cache_ready || !cache.field_cache.isInside(pos)) {
      cache.field_cache_ready = true;
      cache.field_cache       = m_bField.getFieldCell(pos);
    }
    // get the field from the cell
    return std::move(cache.field_cache.getField(pos));
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param[in,out] cache is the propagation cache associated with the track
  ///                      parameters that are being propagated.
  ///
  /// @param[in,out] h     is the desired step size. It can be negative during
  ///                      backwards track propagation, and since we're using an
  ///                      adaptive algorithm, it can also be modified by the
  ///                      stepper class during propagation.
  ///
  double
  step(Cache& cache) const
  {
    // Charge-momentum ratio, in SI units
    const double qop = 1. / units::Nat2SI<units::MOMENTUM>(1. / cache.qop);

    // Runge-Kutta integrator state
    double   h2, half_h;
    Vector3D B_middle, B_last, k2, k3, k4;

    // First Runge-Kutta point (at current position)
    const Vector3D B_first = getField(cache, cache.pos);
    const Vector3D k1      = qop * cache.dir.cross(B_first);

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the local
    // integration error. The results are cached in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory.
    const auto tryRungeKuttaStep = [&](const double h) -> double {
      // Cache the square and half of the step size
      h2     = h * h;
      half_h = h / 2;

      // Second Runge-Kutta point
      const Vector3D pos1 = cache.pos + half_h * cache.dir + h2 / 8 * k1;
      B_middle            = getField(cache, pos1);
      ;
      k2 = qop * (cache.dir + half_h * k1).cross(B_middle);

      // Third Runge-Kutta point
      k3 = qop * (cache.dir + half_h * k2).cross(B_middle);

      // Last Runge-Kutta point
      const Vector3D pos2 = cache.pos + h * cache.dir + h2 / 2 * k3;
      B_last              = getField(cache, pos2);
      k4                  = qop * (cache.dir + h * k3).cross(B_last);

      // Return an estimate of the local integration error
      return h * (k1 - k2 - k3 + k4).template lpNorm<1>();
    };

    // Select and adjust the appropriate Runge-Kutta step size
    double error_estimate = tryRungeKuttaStep(cache.step_size);
    while (error_estimate > 0.0002) {
      cache.step_size *= 0.5;
      error_estimate = tryRungeKuttaStep(cache.step_size);
    }

    // use the adjusted step size
    const double h = cache.step_size;

    // When doing error propagation, update the associated Jacobian matrix
    if (cache.cov_transport) {

      ActsMatrixD<7, 7> D = ActsMatrixD<7, 7>::Zero();
      D(6, 6)             = 1;

      const double conv = units::SI2Nat<units::MOMENTUM>(1);

      auto dFdx = D.block<3, 3>(0, 0);
      auto dFdT = D.block<3, 3>(0, 3);
      auto dFdL = D.block<3, 1>(0, 6);
      auto dGdx = D.block<3, 3>(3, 0);
      auto dGdT = D.block<3, 3>(3, 3);
      auto dGdL = D.block<3, 1>(3, 6);

      ActsMatrixD<3, 3> dk1dT = ActsMatrixD<3, 3>::Zero();
      ActsMatrixD<3, 3> dk2dT = ActsMatrixD<3, 3>::Identity();
      ActsMatrixD<3, 3> dk3dT = ActsMatrixD<3, 3>::Identity();
      ActsMatrixD<3, 3> dk4dT = ActsMatrixD<3, 3>::Identity();

      ActsVectorD<3> dk1dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk2dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk3dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk4dL = ActsVectorD<3>::Zero();

      dk1dL = cache.dir.cross(B_first);
      dk2dL = (cache.dir + half_h * k1).cross(B_middle)
          + qop * half_h * dk1dL.cross(B_middle);
      dk3dL = (cache.dir + half_h * k2).cross(B_middle)
          + qop * half_h * dk2dL.cross(B_middle);
      dk4dL
          = (cache.dir + h * k3).cross(B_last) + qop * h * dk3dL.cross(B_last);

      dk1dT(0, 1) = B_first.z();
      dk1dT(0, 2) = -B_first.y();
      dk1dT(1, 0) = -B_first.z();
      dk1dT(1, 2) = B_first.x();
      dk1dT(2, 0) = B_first.y();
      dk1dT(2, 1) = -B_first.x();
      dk1dT *= qop;

      dk2dT += h / 2 * dk1dT;
      dk2dT = qop * cross(dk2dT, B_middle);

      dk3dT += h / 2 * dk2dT;
      dk3dT = qop * cross(dk3dT, B_middle);

      dk4dT += h * dk3dT;
      dk4dT = qop * cross(dk4dT, B_last);

      dFdx.setIdentity();

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h2 / 6 * (dk1dL + dk2dL + dk3dL);

      dGdx.setZero();

      dGdT.setIdentity();
      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      cache.jacobian = D * cache.jacobian;
    }

    // Update the track parameters according to the equations of motion
    cache.pos += h * cache.dir + h2 / 6 * (k1 + k2 + k3);
    cache.dir += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    cache.dir /= cache.dir.norm();
    cache.derivative.template head<3>()     = cache.dir;
    cache.derivative.template segment<3>(3) = k4;

    // Return the updated step size
    cache.accumulated_path += h;
    return h;
  }

private:
  /// Magnetic field inside of the detector
  BField m_bField;
};

}  // namespace Acts

#endif  // ACTS_EIGEN_STEPPER_HPP
