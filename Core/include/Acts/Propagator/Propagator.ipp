// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/Propagator.hpp"

#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/detail/ActorCalls.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_state_t>
Result<void> Propagator<S, N>::updateNavigationState(
    propagator_state_t& state) const {
  state.position = m_stepper.position(state.stepping);
  const Vector3 tangent = m_stepper.direction(state.stepping);
  state.direction = state.options.direction * tangent;
  state.curvature.setZero();
  const auto curvature = detail::getCurvature(m_stepper, state.stepping);
  if (!curvature.ok()) {
    return curvature.error();
  }
  state.curvature = *curvature;
  return Result<void>::success();
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_state_t>
Result<void> Propagator<S, N>::propagate(propagator_state_t& state) const {
  ACTS_VERBOSE("Entering propagation.");

  state.stage = PropagatorStage::prePropagation;

  // Pre-Propagation: call to the actor, abort condition check
  if (Result<void> preActResult =
          detail::callActorAct(state.options.actor, state.actorState, state,
                               m_stepper, m_navigator, logger());
      !preActResult.ok()) {
    ACTS_DEBUG("Pre-propagation actor call failed: "
               << preActResult.error() << ": "
               << preActResult.error().message());
    return preActResult.error();
  }

  if (detail::callActorCheckAbort(state.options.actor, state.actorState, state,
                                  m_stepper, m_navigator, logger())) {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");

    state.stage = PropagatorStage::postPropagation;

    return detail::callActorAct(state.options.actor, state.actorState, state,
                                m_stepper, m_navigator, logger());
  }

  // Pre-propagation actors may have changed momentum or direction.
  if (auto update = updateNavigationState(state); !update.ok()) {
    return update.error();
  }

  // priming error condition
  bool terminatedNormally = false;

  const auto handleSurfaceReached = [&](const Surface& surface) {
    if constexpr (Concepts::NavigatorHandleSurfaceReachedWithCurvature<N>) {
      m_navigator.handleSurfaceReached(state.navigation, state.position,
                                       state.direction, surface,
                                       state.curvature);
    } else {
      m_navigator.handleSurfaceReached(state.navigation, state.position,
                                       state.direction, surface);
    }
  };

  // Returns none if the actor stops the propagation on a surface it is on
  auto getNextTarget = [&]() -> Result<NavigationTarget> {
    for (unsigned int i = 0; i < state.options.maxTargetSkipping; ++i) {
      const Surface* currentSurface =
          m_navigator.currentSurface(state.navigation);
      NavigationTarget nextTarget = [&] {
        if constexpr (Concepts::NavigatorNextTargetWithCurvature<N>) {
          return m_navigator.nextTarget(state.navigation, state.position,
                                        state.direction, state.curvature);
        } else {
          return m_navigator.nextTarget(state.navigation, state.position,
                                        state.direction);
        }
      }();
      if (nextTarget.isNone()) {
        return NavigationTarget::None();
      }
      IntersectionStatus preStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, nextTarget.surface(), nextTarget.intersectionIndex(),
          state.options.direction, nextTarget.boundaryTolerance(),
          state.options.surfaceTolerance, ConstrainedStep::Type::Navigator,
          logger());
      if (preStepSurfaceStatus == IntersectionStatus::onSurface) {
        // The actors handled the current surface already
        if (&nextTarget.surface() == currentSurface) {
          ACTS_VERBOSE(
              "Pre-step surface status is onSurface, skipping "
              "current surface "
              << currentSurface->geometryId());
          continue;
        }
        // The propagation is on the target already, e.g. at a geometry overlap
        // or at the start. Reach it with a step of zero length.
        ACTS_VERBOSE("Pre-step surface status is onSurface, reaching target "
                     << nextTarget.surface().geometryId() << " without a step");
        state.stage = PropagatorStage::postStep;
        handleSurfaceReached(nextTarget.surface());
        if (Result<void> actResult =
                detail::callActorAct(state.options.actor, state.actorState,
                                     state, m_stepper, m_navigator, logger());
            !actResult.ok()) {
          return actResult.error();
        }
        if (detail::callActorCheckAbort(state.options.actor, state.actorState,
                                        state, m_stepper, m_navigator,
                                        logger())) {
          terminatedNormally = true;
          return NavigationTarget::None();
        }
        // Actors may have changed direction or momentum (and hence curvature).
        if (auto update = updateNavigationState(state); !update.ok()) {
          return update.error();
        }
        state.stage = PropagatorStage::preStep;
        continue;
      }
      if (preStepSurfaceStatus == IntersectionStatus::reachable) {
        return nextTarget;
      }
    }

    ACTS_DEBUG("getNextTarget failed to find a valid target surface after "
               << state.options.maxTargetSkipping << " attempts.");
    return Result<NavigationTarget>::failure(
        PropagatorError::NextTargetLimitReached);
  };

  // Pre-Stepping: target setting
  state.stage = PropagatorStage::preStep;

  Result<NavigationTarget> nextTargetResult = getNextTarget();
  if (!nextTargetResult.ok()) {
    ACTS_DEBUG("Failed to get next target: "
               << nextTargetResult.error() << ": "
               << nextTargetResult.error().message());
    return nextTargetResult.error();
  }
  NavigationTarget nextTarget = *nextTargetResult;

  ACTS_VERBOSE("Starting stepping loop.");

  // Stepping loop
  for (; !terminatedNormally && state.steps < state.options.maxSteps;
       ++state.steps) {
    // Perform a step
    Result<double> res =
        m_stepper.step(state.stepping, state.options.direction,
                       m_navigator.currentVolumeMaterial(state.navigation));
    if (!res.ok()) {
      ACTS_DEBUG("Step failed with " << res.error() << ": "
                                     << res.error().message());
      return res.error();
    }
    // Accumulate the path length
    state.pathLength += *res;
    // Update the position and direction
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    ACTS_VERBOSE("Step with size " << *res << " performed. We are now at "
                                   << state.position.transpose()
                                   << " with direction "
                                   << state.direction.transpose());

    // release actor and aborter constrains after step was performed
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Navigator);
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Actor);

    // Post-stepping: check target status, call actors, check abort conditions
    state.stage = PropagatorStage::postStep;

    if (!nextTarget.isNone()) {
      IntersectionStatus postStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, nextTarget.surface(), nextTarget.intersectionIndex(),
          state.options.direction, nextTarget.boundaryTolerance(),
          state.options.surfaceTolerance, ConstrainedStep::Type::Navigator,
          logger());
      if (postStepSurfaceStatus == IntersectionStatus::onSurface) {
        // A portal transition can create a policy state before actors run.
        if (auto update = updateNavigationState(state); !update.ok()) {
          return update.error();
        }
        handleSurfaceReached(nextTarget.surface());
      }
      if (postStepSurfaceStatus != IntersectionStatus::reachable) {
        nextTarget = NavigationTarget::None();
      }
    }

    Result<void> actResult =
        detail::callActorAct(state.options.actor, state.actorState, state,
                             m_stepper, m_navigator, logger());
    if (!actResult.ok()) {
      ACTS_DEBUG("Actor call failed: " << actResult.error() << ": "
                                       << actResult.error().message());
      return actResult.error();
    }

    if (detail::callActorCheckAbort(state.options.actor, state.actorState,
                                    state, m_stepper, m_navigator, logger())) {
      terminatedNormally = true;
      break;
    }

    // Actors may have changed direction or momentum (and hence curvature).
    if (auto update = updateNavigationState(state); !update.ok()) {
      return update.error();
    }

    // Pre-Stepping: target setting
    state.stage = PropagatorStage::preStep;

    const auto checkTargetValid = [&] {
      if constexpr (Concepts::NavigatorCheckTargetValidWithCurvature<N>) {
        return m_navigator.checkTargetValid(state.navigation, state.position,
                                            state.direction, state.curvature);
      } else {
        return m_navigator.checkTargetValid(state.navigation, state.position,
                                            state.direction);
      }
    };
    if (!nextTarget.isNone() && !checkTargetValid()) {
      ACTS_VERBOSE("Target is not valid anymore.");
      nextTarget = NavigationTarget::None();
    }

    if (nextTarget.isNone()) {
      // navigator step constraint is not valid anymore
      m_stepper.releaseStepSize(state.stepping,
                                ConstrainedStep::Type::Navigator);

      nextTargetResult = getNextTarget();
      if (!nextTargetResult.ok()) {
        ACTS_DEBUG("Failed to get next target: "
                   << nextTargetResult.error() << ": "
                   << nextTargetResult.error().message());
        return nextTargetResult.error();
      }
      if (terminatedNormally) {
        break;
      }
      nextTarget = *nextTargetResult;
    }
  }  // end of stepping loop

  // check if we didn't terminate normally via the actor
  if (!terminatedNormally) {
    ACTS_DEBUG("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << state.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  ACTS_VERBOSE("Stepping loop done.");

  state.stage = PropagatorStage::postPropagation;

  // Post-stepping call to the actor
  if (auto postPropagationResult =
          detail::callActorAct(state.options.actor, state.actorState, state,
                               m_stepper, m_navigator, logger());
      !postPropagationResult.ok()) {
    ACTS_DEBUG("Post-propagation actor call failed: "
               << postPropagationResult.error() << ": "
               << postPropagationResult.error().message());
    return postPropagationResult.error();
  }
  return Result<void>::success();
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_options_t>
auto Propagator<S, N>::propagate(const BoundParameters& start,
                                 const propagator_options_t& options,
                                 bool createFinalParameters) const
    -> Result<ResultType> {
  auto state = makeState(options);
  return propagate(state, start, nullptr, createFinalParameters);
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_options_t>
auto Propagator<S, N>::propagate(const BoundParameters& start,
                                 const Surface& target,
                                 const propagator_options_t& options) const
    -> Result<ResultType> {
  auto state = makeState(options);
  return propagate(state, start, &target, true);
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_state_t>
auto Propagator<S, N>::propagate(propagator_state_t& state,
                                 const BoundParameters& start,
                                 const Surface* target,
                                 bool createFinalParameters) const
    -> Result<ResultType> {
  auto initRes = initialize(state, start, target);
  if (!initRes.ok()) {
    ACTS_DEBUG("Initialization failed: " << initRes.error() << ": "
                                         << initRes.error().message());
    return initRes.error();
  }

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(state, propagationResult, createFinalParameters, target);
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_options_t>
auto Propagator<S, N>::makeState(const propagator_options_t& options) const
    -> State<propagator_options_t> {
  return State<propagator_options_t>{options,
                                     m_stepper.makeState(options.stepping),
                                     m_navigator.makeState(options.navigation)};
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_state_t>
Result<void> Propagator<S, N>::initialize(propagator_state_t& state,
                                          const BoundParameters& start,
                                          const Surface* target) const {
  m_stepper.initialize(state.stepping, start);

  if (auto update = updateNavigationState(state); !update.ok()) {
    return update.error();
  }

  // Navigator initialize state call
  auto navInitRes = m_navigator.initialize(
      state.navigation, {.position = state.position,
                         .direction = state.direction,
                         .curvature = state.curvature,
                         .propagationDirection = state.options.direction,
                         .startSurface = &start.referenceSurface(),
                         .targetSurface = target});
  if (!navInitRes.ok()) {
    ACTS_DEBUG("Navigator initialization failed: "
               << navInitRes.error() << ": " << navInitRes.error().message());
    return navInitRes.error();
  }

  return Result<void>::success();
}

template <StepperConcept S, NavigatorConcept N>
template <typename propagator_state_t>
auto Propagator<S, N>::makeResult(propagator_state_t& state,
                                  Result<void> propagationResult,
                                  bool createFinalParameters,
                                  const Surface* target) const
    -> Result<ResultType> {
  if (!propagationResult.ok()) {
    ACTS_DEBUG("Propagation failed: " << propagationResult.error() << ": "
                                      << propagationResult.error().message());
    return propagationResult.error();
  }

  ResultType result{};
  result.steps = state.steps;
  result.pathLength = state.pathLength;
  result.statistics.stepping = m_stepper.statistics(state.stepping);
  result.statistics.navigation = state.navigation.statistics;

  if (createFinalParameters) {
    if (target == nullptr) {
      target = m_navigator.currentSurface(state.navigation);
    }

    if (target != nullptr) {
      // We are at a surface, so we need to compute the bound state
      const auto jacobian = m_stepper.transportToBound(state.stepping, *target);
      if (!jacobian.ok()) {
        ACTS_DEBUG("Failed to transport to current surface: "
                   << jacobian.error() << ": " << jacobian.error().message());
        return jacobian.error();
      }
      auto parameters = m_stepper.boundParameters(state.stepping, *target);
      if (!parameters.ok()) {
        ACTS_DEBUG("Failed to get bound parameters at current surface: "
                   << parameters.error() << ": "
                   << parameters.error().message());
        return parameters.error();
      }
      result.endParameters = std::move(*parameters);
      if (m_stepper.hasCovariance(state.stepping)) {
        result.transportJacobian = *jacobian;
      }
    } else {
      if (!m_stepper.prepareCurvilinearState(state.stepping)) {
        ACTS_DEBUG("Failed to prepare curvilinear state.");
        return PropagatorError::Failure;
      }
      const auto jacobian = m_stepper.transportToCurvilinear(state.stepping);
      result.endParameters = m_stepper.curvilinearParameters(state.stepping);
      if (m_stepper.hasCovariance(state.stepping)) {
        result.transportJacobian = jacobian;
      }
    }
  }

  return Result<ResultType>::success(std::move(result));
}

template <typename derived_t>
Result<BoundTrackParameters>
detail::BasePropagatorHelper<derived_t>::propagateToSurface(
    const BoundTrackParameters& start, const Surface& target,
    const Options& options) const {
  using DerivedOptions = typename derived_t::template Options<>;

  DerivedOptions derivedOptions(options);

  auto res = static_cast<const derived_t*>(this)->propagate(start, target,
                                                            derivedOptions);

  if (!res.ok()) {
    return res.error();
  }

  // Without errors we can expect a valid endParameters when propagating to a
  // target surface
  assert((*res).endParameters);
  return std::move((*res).endParameters.value());
}

}  // namespace Acts
