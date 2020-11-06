// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/ConstrainedStepControl.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// State for track parameter propagation
///
template <typename stepper_t>
struct StepperState {
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;

  /// Access control to state: Stepper
  friend stepper_t;
  /// Access control to state: StepControl
  friend ConstrainedStepControl<stepper_t>;

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
  template <typename parameters_t>
  explicit StepperState(
      std::reference_wrapper<const GeometryContext> gctx,
      std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
      const parameters_t& par, NavigationDirection ndir = forward,
      double ssize = std::numeric_limits<double>::max(),
      double stolerance = s_onSurfaceTolerance)
      : q(static_cast<int>(par.charge())),
        navDir(ndir),
        stepSize(ndir * std::abs(ssize)),
        tolerance(stolerance),
        geoContext(gctx) {
    pars.template segment<3>(eFreePos0) = par.position(gctx);
    pars.template segment<3>(eFreeDir0) = par.unitDirection();
    pars[eFreeTime] = par.time();
    pars[eFreeQOverP] = par.charge() / par.absoluteMomentum();
    if (par.covariance()) {
      // Get the reference surface for navigation
      const auto& surface = par.referenceSurface();
      // set the covariance transport flag to true and copy
      covTransport = true;
      cov = BoundSymMatrix(*par.covariance());
      jacToGlobal = surface.jacobianLocalToGlobal(gctx, par.parameters());
    }
  }

 protected:
  /// Internal free vector parameters
  FreeVector pars = FreeVector::Zero();

  /// The charge as the free vector can be 1/p or q/p
  int q = 1;

  /// Boolean to indiciate if you need covariance transport
  bool covTransport = false;
  Covariance cov = Covariance::Zero();

  /// Jacobian from local to the global frame
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

  /// Pure transport jacobian part from runge kutta integration
  FreeMatrix jacTransport = FreeMatrix::Identity();

  /// The full jacobian of the transport entire transport
  Jacobian jacobian = Jacobian::Identity();

  /// The propagation derivative
  FreeVector derivative = FreeVector::Zero();

  /// Navigation direction, this is needed for searching
  NavigationDirection navDir;

  /// adaptive step size of the runge-kutta integration
  ConstrainedStep stepSize = std::numeric_limits<double>::max();

  /// The tolerance for the stepping
  double tolerance = s_onSurfaceTolerance;

  // Previous step size for overstep estimation (ignored for SL stepper)
  double previousStepSize = 0.;

  /// accummulated path length state
  double pathAccumulated = 0.;

  // Cache the geometry context of this propagation
  std::reference_wrapper<const GeometryContext> geoContext;
};

}  // namespace Acts
