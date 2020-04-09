// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>
#include <variant>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
  
/// State for track parameter propagation
///
struct StepperState {
  using Jacobian = std::variant<BoundMatrix, FreeToBoundMatrix,
                                BoundToFreeMatrix, FreeMatrix>;
  using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;

  /// Delete the default constructor
  StepperState() = delete;

  /// Constructor from the initial track parameters
  ///
  /// @tparam parameters_t the Type of the track parameters
  ///
  /// @param [in] gctx is the context object for the geometery
  /// @param [in] mctx is the context object for the magnetic field
  /// @param [in] par The track parameters at start
  /// @param [in] ndir is the navigation direction
  /// @param [in] ssize is the (absolute) maximum step size
  /// @param [in] stolerance is the stepping tolerance
  template <typename parameters_t,
            std::enable_if_t<parameters_t::is_local_representation, int> = 0>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max(),
      double stolerance = s_onSurfaceTolerance)
      : pos(par.position()),
        dir(par.momentum().normalized()),
        p(par.momentum().norm()),
        q(par.charge()),
        t(par.time()),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    if (par.covariance()) {
      // Set the covariance transport flag to true
      covTransport = true;
      // Get the covariance
      jacToGlobal = BoundToFreeMatrix::Zero();
      par.referenceSurface().initJacobianToGlobal(gctx, *jacToGlobal, pos, dir,
                                                  par.parameters());
      cov = *par.covariance();
      jacobian.emplace<0>(BoundMatrix::Identity());
    }
  }

  /// Constructor from the initial track parameters
  ///
  /// @tparam parameters_t the Type of the track parameters
  ///
  /// @param [in] gctx is the context object for the geometery
  /// @param [in] mctx is the context object for the magnetic field
  /// @param [in] par The track parameters at start
  /// @param [in] ndir is the navigation direction
  /// @param [in] ssize is the (absolute) maximum step size
  /// @param [in] stolerance is the stepping tolerance
  template <
      typename parameters_t,
      std::enable_if_t<not parameters_t::is_local_representation, int> = 0>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max(),
      double stolerance = s_onSurfaceTolerance)
      : pos(par.position()),
        dir(par.momentum().normalized()),
        p(par.momentum().norm()),
        q(par.charge()),
        t(par.time()),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    if (par.covariance()) {
      // Set the covariance transport flag to true
      covTransport = true;
      // Get the covariance
      cov = *par.covariance();
      jacobian.emplace<3>(FreeMatrix::Identity());
      
      // Set up transformations between angles and directions in jacobian
      jacDirToAngle = directionsToAnglesJacobian(dir);
      jacAngleToDir = anglesToDirectionsJacobian(dir);
    }
  }

  /// Transform from directions to angles in jacobian
  ActsMatrixD<8, 7> jacDirToAngle;
  /// Transform from angles to directions in jacobian
  ActsMatrixD<7, 8> jacAngleToDir;

  /// Jacobian from local to the global frame
  std::optional<BoundToFreeMatrix> jacToGlobal;

  /// Pure transport jacobian part from runge kutta integration
  FreeMatrix jacTransport = FreeMatrix::Identity();

  /// The full jacobian since the first step
  Jacobian jacobian;

  /// The propagation derivative
  FreeVector derivative = FreeVector::Zero();

  /// Boolean to indiciate if you need covariance transport
  bool covTransport = false;
  Covariance cov;

  /// Global particle position
  Vector3D pos = Vector3D(0., 0., 0.);

  /// Momentum direction (normalized)
  Vector3D dir = Vector3D(1., 0., 0.);

  /// Momentum
  double p = 0.;

  /// Save the charge: neutral as default for SL stepper
  double q = 0.;

  /// Propagated time
  double t = 0.;

  /// Navigation direction, this is needed for searching
  NavigationDirection navDir;

  /// accummulated path length state
  double pathAccumulated = 0.;

  /// adaptive step size of the runge-kutta integration
  ConstrainedStep stepSize = std::numeric_limits<double>::max();
  /// Previous step size for overstep estimation (ignored for SL stepper)
  double previousStepSize = 0.;

  /// The tolerance for the stepping
  double tolerance = s_onSurfaceTolerance;

  /// Cache the geometry context of this propagation
  std::reference_wrapper<const GeometryContext> geoContext;
 
 /// TODO: The methods will be move to CovarianceEnginge.cpp & will be called on reinits
/// @brief Constructs a jacobian to transform from (x,y,z,t,Tx,Ty,Tz,q/p) to (x,y,z,t,phi,theta,q/p)
///
/// @param [in] dir Direction vector
///
/// @return The jacobian
  ActsMatrixD<8, 7>
  directionsToAnglesJacobian(const Vector3D dir)
  {
	ActsMatrixD<8, 7> jac = ActsMatrixD<8, 7>::Zero();
	
	const double x = dir(0);  // == cos(phi) * sin(theta)
    const double y = dir(1);  // == sin(phi) * sin(theta)
    const double z = dir(2);  // == cos(theta)
    // can be turned into cosine/sine
    const double cosTheta = z;
    const double sinTheta = sqrt(x * x + y * y);
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = x * invSinTheta;
    const double sinPhi = y * invSinTheta;
    
    jac(0, 0) = 1.;
    jac(1, 1) = 1.;
    jac(2, 2) = 1.;
    jac(3, 3) = 1.;
    jac(7, 6) = 1.;
    
	jac(4, 4) = -sinTheta * sinPhi;
    jac(4, 5) = cosTheta * cosPhi;
    jac(5, 4) = sinTheta * cosPhi;
    jac(5, 5) = cosTheta * sinPhi;
    jac(6, 5) = -sinTheta;
    return jac;
  }
  
  /// @brief Constructs a jacobian to transform from (x,y,z,t,phi,theta,q/p) to (x,y,z,t,Tx,Ty,Tz,q/p)
///
/// @param [in] dir Direction vector
///
/// @return The jacobian
ActsMatrixD<7, 8>
anglesToDirectionsJacobian(const Vector3D dir)
{
ActsMatrixD<7, 8> jacobian = ActsMatrixD<7, 8>::Zero();
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  const double z = dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  jacobian(0, 0) = 1.;
  jacobian(1, 1) = 1.;
  jacobian(2, 2) = 1.;
  jacobian(3, 3) = 1.;
  jacobian(6, 7) = 1.;
  
  jacobian(4, 4) = -sinPhi * invSinTheta;
  jacobian(4, 5) = cosPhi * invSinTheta;
  jacobian(5, 4) = cosPhi * cosTheta;
  jacobian(5, 5) = sinPhi * cosTheta;
  jacobian(5, 6) = -invSinTheta * (1. - cosTheta * cosTheta);

  return jacobian;
}
};
}  // namespace Acts