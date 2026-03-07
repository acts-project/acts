// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/Propagator.hpp"

#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

template <typename S, typename N>
template <typename propagator_state_t>
Result<void> Propagator<S, N>::propagate(propagator_state_t& state) const {
  ACTS_VERBOSE("Entering propagation.");

  state.stage = PropagatorStage::prePropagation;

  // Pre-Propagation: call to the actor list, abort condition check
  state.options.actorList.act(state, m_stepper, m_navigator, logger());

  if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                         logger())) {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");

    state.stage = PropagatorStage::postPropagation;

    return state.options.actorList.act(state, m_stepper, m_navigator, logger());
  }

  auto getNextTarget = [&]() -> Result<NavigationTarget> {
    for (unsigned int i = 0; i < state.options.maxTargetSkipping; ++i) {
      NavigationTarget nextTarget = m_navigator.nextTarget(
          state.navigation, state.position, state.direction);
      if (nextTarget.isNone()) {
        return NavigationTarget::None();
      }
      IntersectionStatus preStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, nextTarget.surface(), nextTarget.intersectionIndex(),
          state.options.direction, nextTarget.boundaryTolerance(),
          state.options.surfaceTolerance, ConstrainedStep::Type::Navigator,
          logger());
      if (preStepSurfaceStatus == IntersectionStatus::onSurface) {
        // This indicates a geometry overlap which is not handled by the
        // navigator, so we skip this target.
        // This can also happen in a well-behaved geometry with external
        // surfaces.
        ACTS_VERBOSE("Pre-step surface status is onSurface, skipping target "
                     << nextTarget.surface().geometryId());
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

  // priming error condition
  bool terminatedNormally = false;

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
  for (; state.steps < state.options.maxSteps; ++state.steps) {
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
        m_navigator.handleSurfaceReached(state.navigation, state.position,
                                         state.direction, nextTarget.surface());
      }
      if (postStepSurfaceStatus != IntersectionStatus::reachable) {
        nextTarget = NavigationTarget::None();
      }
    }

    Result<void> actResult =
        state.options.actorList.act(state, m_stepper, m_navigator, logger());
    if (!actResult.ok()) {
      ACTS_DEBUG("Actor call failed: " << actResult.error() << ": "
                                       << actResult.error().message());
      return actResult.error();
    }

    if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                           logger())) {
      terminatedNormally = true;
      break;
    }

    // Update the position and direction because actors might have changed it
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    // Pre-Stepping: target setting
    state.stage = PropagatorStage::preStep;

    if (!nextTarget.isNone() &&
        !m_navigator.checkTargetValid(state.navigation, state.position,
                                      state.direction)) {
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
      nextTarget = *nextTargetResult;
    }
  }  // end of stepping loop

  // check if we didn't terminate normally via aborters
  if (!terminatedNormally) {
    ACTS_DEBUG("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << state.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  ACTS_VERBOSE("Stepping loop done.");

  state.stage = PropagatorStage::postPropagation;

  // Post-stepping call to the actor list
  if (auto postPropagationResult =
          state.options.actorList.act(state, m_stepper, m_navigator, logger());
      !postPropagationResult.ok()) {
    ACTS_DEBUG("Post-propagation actor call failed: "
               << postPropagationResult.error() << ": "
               << postPropagationResult.error().message());
    return postPropagationResult.error();
  }
  return Result<void>::success();
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Propagator<S, N>::propagate(const parameters_t& start,
                                 const propagator_options_t& options,
                                 bool createFinalParameters) const
    -> Result<ResultType<propagator_options_t>> {
  auto state = makeState<propagator_options_t, path_aborter_t>(options);

  auto initRes =
      initialize<decltype(state), parameters_t, path_aborter_t>(state, start);
  if (!initRes.ok()) {
    ACTS_DEBUG("Initialization failed: " << initRes.error() << ": "
                                         << initRes.error().message());
    return initRes.error();
  }

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, options,
                    createFinalParameters);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Propagator<S, N>::propagate(const parameters_t& start,
                                 const Surface& target,
                                 const propagator_options_t& options) const
    -> Result<ResultType<propagator_options_t>> {
  auto state =
      makeState<propagator_options_t, target_aborter_t, path_aborter_t>(
          target, options);

  auto initRes =
      initialize<decltype(state), parameters_t, path_aborter_t>(state, start);
  if (!initRes.ok()) {
    ACTS_DEBUG("Initialization failed: " << initRes.error() << ": "
                                         << initRes.error().message());
    return initRes.error();
  }

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, options, true,
                    &target);
}

template <typename S, typename N>
template <typename propagator_options_t, typename path_aborter_t>
auto Propagator<S, N>::makeState(const propagator_options_t& options) const {
  // Expand the actor list with a path aborter
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto actorList = options.actorList.append(pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(actorList);

  using OptionsType = decltype(eOptions);
  using StateType = State<OptionsType>;

  StateType state{eOptions, m_stepper.makeState(eOptions.stepping),
                  m_navigator.makeState(eOptions.navigation)};

  return state;
}

template <typename S, typename N>
template <typename propagator_options_t, typename target_aborter_t,
          typename path_aborter_t>
auto Propagator<S, N>::makeState(const Surface& target,
                                 const propagator_options_t& options) const {
  // Expand the actor list with a target and path aborter
  target_aborter_t targetAborter;
  targetAborter.surface = &target;
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto actorList = options.actorList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(actorList);
  eOptions.navigation.targetSurface = &target;

  using OptionsType = decltype(eOptions);
  using StateType = State<OptionsType>;

  StateType state{eOptions, m_stepper.makeState(eOptions.stepping),
                  m_navigator.makeState(eOptions.navigation)};

  return state;
}

template <typename S, typename N>
template <typename propagator_state_t, typename parameters_t,
          typename path_aborter_t>
Result<void> Propagator<S, N>::initialize(propagator_state_t& state,
                                          const parameters_t& start) const {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  m_stepper.initialize(state.stepping, start);

  state.position = m_stepper.position(state.stepping);
  state.direction =
      state.options.direction * m_stepper.direction(state.stepping);

  state.navigation.options.startSurface = &start.referenceSurface();

  // Navigator initialize state call
  auto navInitRes =
      m_navigator.initialize(state.navigation, state.position, state.direction,
                             state.options.direction);
  if (!navInitRes.ok()) {
    ACTS_DEBUG("Navigator initialization failed: "
               << navInitRes.error() << ": " << navInitRes.error().message());
    return navInitRes.error();
  }

  // Apply the loop protection - it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.actorList.template get<path_aborter_t>(),
      false, logger());

  return Result<void>::success();
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Propagator<S, N>::makeResult(
    propagator_state_t state, Result<void> propagationResult,
    const propagator_options_t& /*options*/, bool createFinalParameters,
    const Surface* target) const -> Result<ResultType<propagator_options_t>> {
  // Type of the full propagation result, including output from actors
  using ThisResultType = ResultType<propagator_options_t>;

  if (!propagationResult.ok()) {
    ACTS_DEBUG("Propagation failed: " << propagationResult.error() << ": "
                                      << propagationResult.error().message());
    return propagationResult.error();
  }

  ThisResultType result{};
  moveStateToResult(state, result);

  if (createFinalParameters) {
    if (target == nullptr) {
      target = m_navigator.currentSurface(state.navigation);
    }

    if (target != nullptr) {
      // We are at a surface, so we need to compute the bound state
      if constexpr (!HasMultiStepper) {
        const auto boundState = m_stepper.boundState(state.stepping, *target);
        if (!boundState.ok()) {
          ACTS_DEBUG("Failed to get bound state at current surface: "
                     << boundState.error() << ": "
                     << boundState.error().message());
          return boundState.error();
        }
        result.endParameters = std::get<0>(*boundState);
        if (state.stepping.covTransport) {
          result.transportJacobian = std::get<1>(*boundState);
        }
      } else {
        const auto multiBoundState =
            m_stepper.multiBoundState(state.stepping, *target);
        if (!multiBoundState.ok()) {
          ACTS_DEBUG("Failed to get multi bound state at current surface: "
                     << multiBoundState.error() << ": "
                     << multiBoundState.error().message());
          return multiBoundState.error();
        }
        result.endParameters =
            std::get<0>(*multiBoundState)
                .merge(state.stepping.options.componentMergeMethod);
        result.endParametersMultiComponent = std::get<0>(*multiBoundState);
        if (state.stepping.covTransport) {
          result.transportJacobian = std::get<1>(*multiBoundState);
        }
      }
    } else {
      if (!m_stepper.prepareCurvilinearState(state.stepping)) {
        ACTS_DEBUG("Failed to prepare curvilinear state.");
        return PropagatorError::Failure;
      }
      if constexpr (!HasMultiStepper) {
        const auto curvState = m_stepper.curvilinearState(state.stepping);
        result.endParameters = std::get<BoundTrackParameters>(curvState);
        if (state.stepping.covTransport) {
          result.transportJacobian = std::get<1>(curvState);
        }
      } else {
        const auto multiCurvState =
            m_stepper.multiCurvilinearState(state.stepping);
        result.endParameters =
            std::get<0>(multiCurvState)
                .merge(state.stepping.options.componentMergeMethod);
        result.endParametersMultiComponent = std::get<0>(multiCurvState);
        if (state.stepping.covTransport) {
          result.transportJacobian = std::get<1>(multiCurvState);
        }
      }
    }
  }

  return Result<ThisResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_result_t>
void Propagator<S, N>::moveStateToResult(propagator_state_t& state,
                                         propagator_result_t& result) const {
  result.tuple() = std::move(state.tuple());

  result.steps = state.steps;
  result.pathLength = state.pathLength;

  result.statistics.stepping = state.stepping.statistics;
  result.statistics.navigation = state.navigation.statistics;
}

template <typename derived_t>
Result<BoundTrackParameters>
detail::BasePropagatorHelper<derived_t>::propagateToSurface(
    const BoundTrackParameters& start, const Surface& target,
    const Options& options) const {
  using DerivedOptions = typename derived_t::template Options<>;
  using DerivedResult = typename derived_t::template ResultType<DerivedOptions>;

  DerivedOptions derivedOptions(options);

  // dummy initialization
  Result<DerivedResult> res =
      Result<DerivedResult>::failure(PropagatorError::Failure);

  // Due to the geometry of the perigee surface the overstepping tolerance
  // is sometimes not met.
  if (target.type() == Surface::SurfaceType::Perigee) {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   ForcedSurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  } else {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   SurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  }

  if (!res.ok()) {
    return res.error();
  }

  // Without errors we can expect a valid endParameters when propagating to a
  // target surface
  assert(res->endParameters &&
         "Expected valid end parameters when propagating to a target surface.");
  return std::move(res->endParameters.value());
}

}  // namespace Acts
