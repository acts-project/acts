// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/detail/NavigationStepGuard.hpp"

namespace Acts::detail {
/// Navigator-specific part of the optional propagation-step adapter. The
/// actual state lives with navigation; geometry is supplied by its policies.
template <>
struct NavigationStepGuard<Navigator> {
  template <typename Stepper, typename State>
  static void beforeStep(const Stepper& stepper, State& state,
                         NavigationTarget& nextTarget) {
    // The certificate follows one trajectory and needs its individual step
    // constraints. A multi-stepper's aggregate state cannot provide this.
    if constexpr (Concepts::SingleStepper<Stepper>) {
      if (!state.navigation.options.experimentalCylindricalPortalGuard) {
        return;
      }
      auto& portalGuard = state.navigation.cylinderPortalGuard;
      portalGuard.positiveStep = false;
      const bool runPortalGuard = portalGuard.needsCheck(
          state.options.geoContext, state.navigation, state.position,
          state.direction, state.curvature, nextTarget,
          stepper.stepSize(state.stepping),
          state.options.direction * state.pathLength);
      if (runPortalGuard) {
        portalGuard.beforeStep(state.options.geoContext, state.navigation,
                               state.position, state.direction,
                               stepper.stepSize(state.stepping).value());
        if (portalGuard.pendingPortal != nullptr) {
          stepper.updateStepSize(state.stepping, portalGuard.pendingFlight,
                                 ConstrainedStep::Type::Actor);
        }
        const auto radialTarget = portalGuard.radialAlternative(
            state.options.geoContext, state.navigation, state.position,
            state.direction, nextTarget,
            stepper.stepSize(state.stepping).value(),
            stepper.stepSize(state.stepping)
                .value(ConstrainedStep::Type::Navigator));
        if (!radialTarget.isNone()) {
          // Select the actual native cylinder root before stepping. Its
          // subsequent negative convergence steps and normal actor arrival
          // remain untouched.
          state.navigation.pendingTarget.reset();
          nextTarget = radialTarget;
          state.navigation.stream.currentCandidate() = radialTarget;
          portalGuard.radialOverride = true;
          stepper.releaseStepSize(state.stepping,
                                  ConstrainedStep::Type::Navigator);
          stepper.updateStepSize(state.stepping, radialTarget.pathLength(),
                                 ConstrainedStep::Type::Navigator);
        }
        portalGuard.positiveStep =
            stepper.stepSize(state.stepping).value() > 0.;
        const double radialCap = portalGuard.radialStepLimit(
            state.navigation, state.position, state.direction, state.curvature,
            nextTarget, stepper.stepSize(state.stepping).value());
        if (portalGuard.positiveStep) {
          stepper.updateStepSize(state.stepping, radialCap,
                                 ConstrainedStep::Type::Actor);
        }
      }
    }
  }

  template <typename State>
  static void afterStep(State& state, NavigationTarget& nextTarget,
                        double performedStep) {
    if (!state.navigation.options.experimentalCylindricalPortalGuard) {
      return;
    }
    auto& portalGuard = state.navigation.cylinderPortalGuard;
    // A radial face can be coincident within tolerance with the selected disc.
    // Stage the actual outward radial exit before ordinary post-step handling;
    // do not silently skip it as an entry surface or switch volume in
    // Navigator.
    const Portal* physicalRadialPortal = portalGuard.physicalRadialExit(
        state.options.geoContext, state.navigation, state.position,
        state.direction, portalGuard.positiveStep);
    if (physicalRadialPortal != nullptr &&
        (!nextTarget.isPortalTarget() ||
         &nextTarget.portal() != physicalRadialPortal)) {
      const auto intersections = physicalRadialPortal->surface().intersect(
          state.options.geoContext, state.position, state.direction,
          BoundaryTolerance::None(), state.options.surfaceTolerance);
      const auto [intersection, index] = intersections.closestWithIndex();
      if (intersection.status() == IntersectionStatus::onSurface) {
        state.navigation.pendingTarget.reset();
        nextTarget =
            NavigationTarget(intersection, index, *physicalRadialPortal,
                             BoundaryTolerance::None());
        state.navigation.stream.currentCandidate() = nextTarget;
      }
    }

    // Promote only an actual finite owner-face crossing to a Portal arrival.
    // The ordinary post-step handler and actor sequence below delivers it;
    // there is no direct Navigator-only transition or skipped portal material.
    if (physicalRadialPortal == nullptr &&
        portalGuard.isPhysicalExit(state.options.geoContext, state.navigation,
                                   state.position, state.direction,
                                   performedStep)) {
      // An additional surface can have displaced a geometry candidate. A
      // portal arrival changes volume, so that saved target is obsolete.
      state.navigation.pendingTarget.reset();
      nextTarget = NavigationTarget(
          Intersection3D(state.position, 0., IntersectionStatus::onSurface), 0,
          *portalGuard.pendingPortal, BoundaryTolerance::None());
      state.navigation.stream.currentCandidate() = nextTarget;
    }
  }
  template <typename State>
  static void targetStatus(State& state,
                           IntersectionStatus postStepSurfaceStatus) {
    if (!state.navigation.options.experimentalCylindricalPortalGuard) {
      return;
    }
    auto& portalGuard = state.navigation.cylinderPortalGuard;
    if (portalGuard.radialOverride &&
        postStepSurfaceStatus == IntersectionStatus::unreachable) {
      // A local ray root need not remain reachable along a curved flight.
      // Ask the original policy for a fresh pool after such a rejection;
      // do not consume the overwritten old candidate as if it were reached.
      state.navigation.stream.reset();
    }
    portalGuard.radialOverride = false;
  }
};
}  // namespace Acts::detail
