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

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// State for track parameter propagation
///
struct StepperState {
  using cstep = detail::ConstrainedStep;

  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;

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
  template <typename parameters_t>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max())
      : pos(par.position()),
        dir(par.momentum().normalized()),
        p(par.momentum().norm()),
        q(par.charge()),
        t0(par.time()),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    if (par.covariance()) {  // TODO: constructors might be combined but a
                             // covariance setter is then templated
      // Set the covariance transport flag to true
      covTransport = true;
      // Get the covariance
      par.referenceSurface().initJacobianToGlobal(gctx, jacToGlobal, pos, dir,
                                                  par.parameters());
      cov = *par.covariance();
    }
  }

  /// Jacobian from local to the global frame
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

  /// Pure transport jacobian part from runge kutta integration
  FreeMatrix jacTransport = FreeMatrix::Identity();

  /// The full jacobian since the first step
  Jacobian jacobian = Jacobian::Identity();

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

  /// @note The time is split into a starting and a propagated time to avoid
  /// machine precision related errors
  /// Starting time
  const double t0;
  /// Propagated time
  double dt = 0.;

  /// Navigation direction, this is needed for searching
  NavigationDirection navDir;

  /// accummulated path length state
  double pathAccumulated = 0.;

  /// adaptive step size of the runge-kutta integration
  cstep stepSize = std::numeric_limits<double>::max();

  /// Previous step size for overstep estimation (ignored for SL stepper)
  double previousStepSize = 0.;

  /// The tolerance for the stepping
  double tolerance = s_onSurfaceTolerance;

  /// Cache the geometry context of this propagation
  std::reference_wrapper<const GeometryContext> geoContext;
};
}  // namespace Acts