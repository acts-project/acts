// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
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
  /// Cache for track parameter propagation
  ///
  struct Cache
  {
    /// Constructor from the initial track parameters
    /// @param [in] par The track parameters at start
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
      // Init the jacobian matrix if needed
      if (par.covariance()) {
        cov_transport = true;
        cov           = ActsSymMatrixD<5>(*par.covariance());
        // Get the reference surface
        const auto& surface = par.referenceSurface();
        surface.initJacobianToGlobal(jac_to_global, pos, dir, par.parameters());
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

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the cache
    ///
    /// @param reinitialize is a flag to steer whether the
    ///        cache should be reinitialized at the new
    ///        position
    ///
    /// @return the full transport jacobian for an eveutal observer
    const ActsMatrixD<5, 5>
    apply_cov_transport(bool reinitialize = false)
    {
      // Optimized trigonometry on the propagation direction
      const double x = dir(0);  // == cos(phi) * sin(theta)
      const double y = dir(1);  // == sin(phi) * sin(theta)
      const double z = dir(2);  // == cos(theta)
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
      // Apply the transport from the steps on the jacobian
      jac_to_global = jac_transport * jac_to_global;
      // Transport the covariance
      ActsRowVectorD<3>       norm_vec(dir);
      const ActsRowVectorD<5> sfactors
          = norm_vec * jac_to_global.topLeftCorner<3, 5>();
      // The full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jac_full
          = jac_to_curv * (jac_to_global - derivative * sfactors);
      // Apply the actual covariance transport
      cov = (jac_full * cov * jac_full.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jac_to_global = ActsMatrixD<7, 5>::Zero();
        jac_transport = ActsMatrixD<7, 7>::Identity();
        // fill the jacobian to global for next transport
        jac_to_global(0, eLOC_0) = -sin_phi;
        jac_to_global(0, eLOC_1) = -cos_phi * cos_theta;
        jac_to_global(1, eLOC_0) = cos_phi;
        jac_to_global(1, eLOC_1) = -sin_phi * cos_theta;
        jac_to_global(2, eLOC_1) = sin_theta;
        jac_to_global(3, ePHI)   = -sin_theta * sin_phi;
        jac_to_global(3, eTHETA) = cos_theta * cos_phi;
        jac_to_global(4, ePHI)   = sin_theta * cos_phi;
        jac_to_global(4, eTHETA) = cos_theta * sin_phi;
        jac_to_global(5, eTHETA) = -sin_theta;
        jac_to_global(6, eQOP)   = 1;
      }
      // return the full transport jacobian (for an observer)
      return jac_full;
    }

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the cache
    ///
    /// @tparam S the Surfac type
    ///
    /// @param surface is the surface to which the covariance is
    ///        forwarded to
    /// @param reinitialize is a flag to steer whether the
    ///        cache should be reinitialized at the new
    ///        position
    /// @note no check is done if the position is actually on the surface
    ///
    /// @return the full transport jacobian for an eveutal observer
    template <typename S>
    const ActsMatrixD<5, 5>
    apply_cov_transport(const S& surface, bool reinitialize = false)
    {
      // Initialize the transport final frame jacobian
      ActsMatrixD<5, 7> jac_to_local = ActsMatrixD<5, 7>::Zero();
      // initalize the jacobian to local, returns the transposed ref frame
      auto rframeT = surface.initJacobianToLocal(jac_to_local, pos, dir);
      // Update the jacobian with the transport from the steps
      jac_to_global = jac_transport * jac_to_global;
      // calculate the form factors for the derivatives
      const ActsRowVectorD<5> s_vec
          = surface.derivativeFactors(pos, dir, rframeT, jac_to_global);
      // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jac_full
          = jac_to_local * (jac_to_global - derivative * s_vec);
      // Apply the actual covariance transport
      cov = (jac_full * cov * jac_full.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jac_to_global = ActsMatrixD<7, 5>::Zero();
        jac_transport = ActsMatrixD<7, 7>::Identity();
        // reset the derivative
        derivative = ActsVectorD<7>::Zero();
        // fill the jacobian to global for next transport
        Vector2D loc{0., 0.};
        surface.globalToLocal(pos, dir, loc);
        ActsVectorD<5> pars;
        pars << loc[eLOC_0], loc[eLOC_1], dir.phi(), dir.theta(), qop;
        surface.initJacobianToGlobal(jac_to_global, pos, dir, pars);
      }
      // return the full transport jacobian (for an observer)
      return jac_full;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);

    /// Charge-momentum ratio, in natural units
    double qop = 1;

    /// Jacobian from local to the global frame
    ActsMatrixD<7, 5> jac_to_global = ActsMatrixD<7, 5>::Zero();

    /// Pure transport jacobian part from runge kutta integration
    ActsMatrixD<7, 7> jac_transport = ActsMatrixD<7, 7>::Identity();

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
  /// - when propagating to a surface we usually return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor requires knowledge of the detector's magnetic field
  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  /// Convert the propagation cache (global) to curvilinear parameters
  /// @param cache The stepper cache
  /// @param reinitialize is a flag to (optionally) reinitialse the cache
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(Cache& cache, bool reinitialize = false)
  {
    double                                   charge = cache.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov_ptr = nullptr;
    // only do the covariance transport if needed
    if (cache.cov_transport) {
      // transport the covariance forward
      cache.apply_cov_transport(reinitialize);
      cov_ptr = std::make_unique<const ActsMatrixD<5, 5>>(cache.cov);
    }
    // return the parameters
    return CurvilinearParameters(
        std::move(cov_ptr), cache.pos, cache.dir / std::abs(cache.qop), charge);
  }

  /// Convert the propagation cache to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] cache Propagation cache used
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename S>
  static BoundParameters
  convert(Cache& cache, const S& surface, bool reinitialize = false)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> cov_ptr = nullptr;
    // Perform error propagation if an initial covariance matrix was provided
    if (cache.cov_transport) {
      // transport the covariance forward
      cache.apply_cov_transport(surface, reinitialize);
      cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cache.cov);
    }
    double charge = cache.qop > 0. ? 1. : -1.;
    // return the bound parameters
    return BoundParameters(std::move(cov_ptr),
                           cache.pos,
                           cache.dir / std::abs(cache.qop),
                           charge,
                           surface);
  }

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
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
    // allowing integration to continue once the error is deemed satisfactory
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
    // @todo remove magic numbers and implement better step estimation
    double error_estimate = tryRungeKuttaStep(cache.step_size);
    while (error_estimate > 0.0002) {
      cache.step_size *= 0.5;
      error_estimate = tryRungeKuttaStep(cache.step_size);
    }

    // use the adjusted step size
    const double h = cache.step_size;

    // When doing error propagation, update the associated Jacobian matrix
    if (cache.cov_transport) {
      // The step transport matrix in global coordinates
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

      // for moment, only update the transport part
      cache.jac_transport = D * cache.jac_transport;
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
