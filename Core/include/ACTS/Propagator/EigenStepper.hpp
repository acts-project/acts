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

public:
  /// enumeration definition of transport parameters
  /// this aims to make the jacobian be read more easily
  enum tPAR { tX0 = 0, tX1 = 1, tX2 = 2, tD0 = 3, tD1 = 4, tD2 = 5, tQOP = 6 };

  /// Cache for track parameter propagation
  ///
  /// it is exposed to public for use of the expert-only
  /// propagate_with_cache method of the propagator
  struct Cache
  {
    /// Constructor from the initial track parameters
    /// @param [in] par The track parameters at start
    ///
    /// @note the covariance matrix is copied when needed
    template <typename T>
    explicit Cache(const T& par)
      : cache_ready(false)
      , pos(par.position())
      , dir(par.momentum().normalized())
      , qop(par.charge() / par.momentum().norm())
      , cov_transport(false)
      , accumulated_path(0.)
      , step_size(std::numeric_limits<double>::max())
    {
      init_jacobian(par);
      // the cache is ready now
      cache_ready = true;
    }

    /// The cache update for optimal performance
    /// @param [in] par The new track parameters at start
    ///
    /// @todo check to identify an reuse of start/cache
    template <typename T>
    void
    update(const T& par)
    {
      if (!cache_ready) {
        pos           = par.position();
        dir           = par.momentum().normalized();
        qop           = (par.charge() / par.momentum().norm());
        cov_transport = false;
        init_jacobian(par);
        cache_ready = true;
      }
    }

    /// The jacobiann initialisation
    /// @param [in] par The (new) track parameters at start
    template <typename T>
    void
    init_jacobian(const T& par)
    {
      if (par.covariance() && !cache_ready) {
        cov_transport = true;
        cov           = ActsSymMatrixD<5>(*par.covariance());
        // The trigonometry required to convert the direction to spherical
        // coordinates and then compute the sines and cosines again can be
        // surprisingly expensive from a performance point of view.
        //
        // Here, we can avoid it because the direction is by definition a unit
        // vector, with the following coordinate conversions...
        const double x = dir(0);  // == cos(phi) * sin(theta)
        const double y = dir(1);  // == sin(phi) * sin(theta)
        const double z = dir(2);  // == cos(theta)

        // ...which we can invert to directly get the sines and cosines:
        const double cos_theta     = z;
        const double sin_theta     = sqrt(x * x + y * y);
        const double inv_sin_theta = 1. / sin_theta;
        const double cos_phi       = x * inv_sin_theta;
        const double sin_phi       = y * inv_sin_theta;

        // get the reference surface
        const auto& surface = par.referenceSurface();
        // the error is (usually) given on a planar-type surface
        // called the reference frame
        const auto& transform = par.referenceFrame();
        // the local error components - it's given by the reference
        // frame for all surfaces but a Disc surface
        if (surface.type() != Surface::Disc) {
          jacobian(tX0, eLOC_0) = transform(0, 0);
          jacobian(tX0, eLOC_1) = transform(0, 1);
          jacobian(tX1, eLOC_0) = transform(1, 0);
          jacobian(tX1, eLOC_1) = transform(1, 1);
          jacobian(tX2, eLOC_0) = transform(2, 0);
          jacobian(tX2, eLOC_1) = transform(2, 1);
        } else {
          // special polar coordinates for the Disc
          const auto& parameters = par.parameters();
          double      lrad       = parameters[eLOC_0];
          double      lphi       = parameters[eLOC_1];
          double      lcos_phi   = cos(lphi);
          double      lsin_phi   = sin(lphi);
          // fill the jacobian entries in block form
          jacobian.template block<3, 1>(tX0, eLOC_0)
              = lcos_phi * transform.template block<3, 1>(0, 0)
              + lsin_phi * transform.template block<3, 1>(0, 1);
          jacobian.template block<3, 1>(tX0, eLOC_1)
              = lrad * (lcos_phi * transform.template block<3, 1>(0, 1)
                        - lsin_phi * transform.template block<3, 1>(0, 0));
        }
        // the momentum components
        jacobian(tD0, ePHI)   = (-sin_theta) * sin_phi;
        jacobian(tD0, eTHETA) = cos_theta * cos_phi;
        jacobian(tD1, ePHI)   = sin_theta * cos_phi;
        jacobian(tD1, eTHETA) = cos_theta * sin_phi;
        jacobian(tD2, eTHETA) = (-sin_theta);
        jacobian(tQOP, eQOP)  = 1;
        // special components for line and straw
        if (surface.type() == Surface::Perigee
            || surface.type() == Surface::Straw) {
          // get the parameter vector
          const auto& pv = par.parameters();
          // the projection of direction onto ref frame normal
          double ipdn = 1. / dir.dot(transform.col(2));
          // build the cross product of d(D)/d(ePHI) components with y axis
          auto dDPhiY = transform.template block<3, 1>(0, 1).template cross(
              jacobian.template block<3, 1>(tD0, ePHI));
          // and the same for the d(D)/d(eTheta) components
          auto dDThetaY = transform.template block<3, 1>(0, 1).template cross(
              jacobian.template block<3, 1>(tD0, eTHETA));
          // and correct for the x axis components
          dDPhiY -= transform.template block<3, 1>(0, 0)
              * (transform.template block<3, 1>(0, 0).template dot(dDPhiY));
          dDThetaY -= transform.template block<3, 1>(0, 0)
              * (transform.template block<3, 1>(0, 0).template dot(dDThetaY));
          // set the jacobian components
          jacobian.template block<3, 1>(tX0, ePHI) = dDPhiY * pv[eLOC_0] * ipdn;
          jacobian.template block<3, 1>(tX0, eTHETA)
              = dDThetaY * pv[eLOC_0] * ipdn;
        }
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

    /// Cache ready flag, avoids double inialization
    bool cache_ready = false;

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
    /// It caches the current magnetic field cell and stays interpolates within
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

  /// Convert the propagation cache (global) to curvilinear parameters
  /// @param cache The stepper cache
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(Cache& cache)
  {
    double                                   charge = cache.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov    = nullptr;

    if (cache.cov_transport) {
      // Optimized trigonometry on the propagation direction, see documentation
      // of Cache::init_jacobian() for a longer mathematical discussion.
      const double x = cache.dir(0);  // == cos(phi) * sin(theta)
      const double y = cache.dir(1);  // == sin(phi) * sin(theta)
      const double z = cache.dir(2);  // == cos(theta)
      // can be turned into cosine/sine
      const double cos_theta     = z;
      const double sin_theta     = sqrt(x * x + y * y);
      const double inv_sin_theta = 1. / sin_theta;
      const double cos_phi       = x * inv_sin_theta;
      const double sin_phi       = y * inv_sin_theta;
      // prepare the jacobian to curvilinear
      ActsMatrixD<5, 7> jac_to_curv = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cos_theta) < s_curvilinearProjTolerance) {
        // We normally operate in curvilinear coordinates defined as follows
        jac_to_curv(0, 0) = -sin_phi;
        jac_to_curv(0, 1) = cos_phi;
        jac_to_curv(1, 0) = -cos_phi * cos_theta;
        jac_to_curv(1, 1) = -sin_phi * cos_theta;
        jac_to_curv(1, 2) = sin_theta;
      } else {
        // Under grazing incidence to z, the above coordinate system definition
        // becomes numerically unstable, and we need to switch to another one
        const double c     = sqrt(y * y + z * z);
        const double inv_c = 1. / c;
        jac_to_curv(0, 1) = -z * inv_c;
        jac_to_curv(0, 2) = y * inv_c;
        jac_to_curv(1, 0) = c;
        jac_to_curv(1, 1) = -x * y * inv_c;
        jac_to_curv(1, 2) = -x * z * inv_c;
      }
      // Directional and momentum parameters for curvilinear
      jac_to_curv(2, 3) = -sin_phi * inv_sin_theta;
      jac_to_curv(2, 4) = cos_phi * inv_sin_theta;
      jac_to_curv(3, 5) = -inv_sin_theta;
      jac_to_curv(4, 6) = 1;

      // transport the covariance
      ActsRowVectorD<3> norm_vec(cache.dir);
      auto cov_at_frame = transport_covariance(cache, jac_to_curv, norm_vec);
      // create the new covariance matrix: curvilinear == measurement frame
      cov = std::make_unique<const ActsSymMatrixD<5>>(std::move(cov_at_frame));
    }
    // this invalidates the cache
    cache.cache_ready = false;
    // return the parameters
    return CurvilinearParameters(
        std::move(cov), cache.pos, cache.dir / std::abs(cache.qop), charge);
  }

  /// Convert the propagation cache to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] cache Propagation cache used
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename S>
  static BoundParameters
  convert(Cache& cache, const S& surface)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> cov = nullptr;

    // Perform error propagation if an initial covariance matrix was provided
    if (cache.cov_transport) {
      // Initialize the transport final frame jacobian
      ActsMatrixD<5, 7> jac_to_local = ActsMatrixD<5, 7>::Zero();
      // Optimized trigonometry on the propagation direction, see documentation
      // of Cache::init_jacobian() for a longer mathematical discussion.
      const double x = cache.dir(0);  // == cos(phi) * sin(theta)
      const double y = cache.dir(1);  // == sin(phi) * sin(theta)
      //
      const double inv_sin_theta_2        = 1. / (x * x + y * y);
      const double cos_phi_over_sin_theta = x * inv_sin_theta_2;
      const double sin_phi_over_sin_theta = y * inv_sin_theta_2;
      const double inv_sin_theta          = sqrt(inv_sin_theta_2);
      // The measurement frame of the surface
      // @todo deal with the disc surface
      RotationMatrix3D rframeT
          = surface.referenceFrame(cache.pos, cache.dir).transpose();
      jac_to_local.block<2, 3>(0, 0) = rframeT.template block<2, 3>(0, 0);
      // Directional and momentum elements for reference frame surface
      jac_to_local(ePHI, 3)   = -sin_phi_over_sin_theta;
      jac_to_local(ePHI, 4)   = cos_phi_over_sin_theta;
      jac_to_local(eTHETA, 5) = -inv_sin_theta;
      jac_to_local(eQOP, 6)   = 1;
      // Create the normal and scale it with the projection onto the direction
      ActsRowVectorD<3> norm_vec = rframeT.template block<1, 3>(2, 0);
      norm_vec /= (norm_vec * cache.dir);
      // calculate the jacobian at the measurement frame
      auto cov_at_frame = transport_covariance(cache, jac_to_local, norm_vec);
      cov = std::make_unique<const ActsSymMatrixD<5>>(std::move(cov_at_frame));
    }
    double charge = cache.qop > 0. ? 1. : -1.;
    // this invalidates the cache
    cache.cache_ready = false;
    // return the bound parameters
    return BoundParameters(std::move(cov),
                           cache.pos,
                           cache.dir / std::abs(cache.qop),
                           charge,
                           surface);
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
  ///                      the cache contains the desired step size.
  ///                      It can be negative during ackwards track propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      also
  ///                      be modified by the stepper class during propagation.
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
      k2                  = qop * (cache.dir + half_h * k1).cross(B_middle);

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
    // @todo remove magic numbers
    double error_estimate = tryRungeKuttaStep(cache.step_size);
    while (error_estimate > 0.0002) {
      cache.step_size *= 0.5;
      error_estimate = tryRungeKuttaStep(cache.step_size);
    }

    // use the adjusted step size
    const double h = cache.step_size;

    // When doing error propagation, update the associated Jacobian matrix
    if (cache.cov_transport) {

      ActsMatrixD<7, 7> D = ActsMatrixD<7, 7>::Identity();
      const double conv = units::SI2Nat<units::MOMENTUM>(1);

      // This sets the reference to the sub matrices
      // dFdx is already initialised as (3x3) idendity
      auto dFdT = D.block<3, 3>(0, 3);
      auto dFdL = D.block<3, 1>(0, 6);
      // dGdx is already initialised as (3x3) zero
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

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h2 / 6 * (dk1dL + dk2dL + dk3dL);

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
  /// @brief Convert the covariance matrix at the local frame
  ///
  /// @param cache The propagation cache
  /// @parm jac_local The Jacobian to Local
  /// @param norm_vec The normal vector
  ///
  /// @return a 5x5 covariance matrix after transport
  static ActsSymMatrixD<5>
  transport_covariance(const Cache& cache,
                       const ActsMatrixD<5, 7>& jac_local,
                       const ActsRowVectorD<3>& norm_vec)
  {
    const ActsRowVectorD<5> scale_factors
        = norm_vec * cache.jacobian.template topLeftCorner<3, 5>();
    // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
    const ActsMatrixD<5, 5> jac
        = jac_local * (cache.jacobian - cache.derivative * scale_factors);
    // return the transported and local covariance matrix
    return ActsSymMatrixD<5>(jac * cache.cov * jac.transpose());
  }

  /// Magnetic field inside of the detector
  BField m_bField;
};

}  // namespace Acts

#endif  // ACTS_EIGEN_STEPPER_HPP
