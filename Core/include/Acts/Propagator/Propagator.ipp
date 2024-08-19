// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <type_traits>

template <typename S, typename N>
template <typename propagator_state_t>
auto Acts::Propagator<S, N>::propagate(propagator_state_t& state) const
    -> Result<void> {
  // Pre-stepping call to the navigator and action list
  ACTS_VERBOSE("Entering propagation.");

  state.stage = PropagatorStage::prePropagation;

  // Pre-Stepping call to the action list
  state.options.actionList(state, m_stepper, m_navigator, logger());
  // assume negative outcome, only set to true later if we actually have
  // a positive outcome.

  // start at true, if we don't begin the stepping loop we're fine.
  bool terminatedNormally = true;

  // Pre-Stepping: abort condition check
  if (state.options.abortList(state, m_stepper, m_navigator, logger())) {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");
  } else {
    // Stepping loop
    ACTS_VERBOSE("Starting stepping loop.");

    // priming error condition
    terminatedNormally = false;

    SurfaceIntersection nextTargetIntersection = SurfaceIntersection::invalid();

    // Propagation loop : stepping
    for (; state.steps < state.options.maxSteps; ++state.steps) {
      // Pre-Stepping: target setting
      state.stage = PropagatorStage::preStep;

      if (!nextTargetIntersection.isValid()) {
        for (int i = 0; i < 3; ++i) {
          nextTargetIntersection = m_navigator.estimateNextTarget(
              state.navigation, state.position,
              state.options.direction * state.direction);
          if (!nextTargetIntersection.isValid()) {
            break;
          }
          IntersectionStatus preStepSurfaceStatus =
              m_stepper.updateSurfaceStatus(
                  state.stepping, *nextTargetIntersection.object(),
                  nextTargetIntersection.index(), state.options.direction,
                  BoundaryTolerance::None(), s_onSurfaceTolerance, logger());
          if (preStepSurfaceStatus >= Acts::IntersectionStatus::reachable) {
            break;
          }
          m_navigator.registerSurfaceStatus(
              state.navigation, state.position,
              state.options.direction * state.direction,
              *nextTargetIntersection.object(), preStepSurfaceStatus);
        }
      }

      // Perform a propagation step - it takes the propagation state
      Result<double> res = m_stepper.step(state, m_navigator);
      if (!res.ok()) {
        ACTS_ERROR("Step failed with " << res.error() << ": "
                                       << res.error().message());
        // pass error to caller
        return res.error();
      }
      // Accumulate the path length
      state.pathLength += *res;
      // Update the position and direction
      state.position = m_stepper.position(state.stepping);
      state.direction = m_stepper.direction(state.stepping);

      ACTS_VERBOSE("Step with size = " << *res << " performed. We are now at: "
                                       << state.position.transpose()
                                       << " with direction: "
                                       << state.direction.transpose());

      // release actor and aborter constrains after step was performed
      m_stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);
      m_stepper.releaseStepSize(state.stepping, ConstrainedStep::aborter);

      // Post-stepping:
      // navigator - action list - aborter list
      state.stage = PropagatorStage::postStep;

      if (nextTargetIntersection.isValid()) {
        IntersectionStatus postStepSurfaceStatus =
            m_stepper.updateSurfaceStatus(
                state.stepping, *nextTargetIntersection.object(),
                nextTargetIntersection.index(), state.options.direction,
                BoundaryTolerance::None(), s_onSurfaceTolerance, logger());
        m_navigator.registerSurfaceStatus(
            state.navigation, state.position,
            state.options.direction * state.direction,
            *nextTargetIntersection.object(), postStepSurfaceStatus);
        if (postStepSurfaceStatus == IntersectionStatus::onSurface) {
          nextTargetIntersection = SurfaceIntersection::invalid();
        }
      }

      state.options.actionList(state, m_stepper, m_navigator, logger());

      if (state.options.abortList(state, m_stepper, m_navigator, logger())) {
        terminatedNormally = true;
        break;
      }
    }
  }

  state.stage = PropagatorStage::postPropagation;

  // if we didn't terminate normally (via aborters) set navigation break.
  // this will trigger error output in the lines below
  if (!terminatedNormally) {
    m_navigator.navigationBreak(state.navigation, true);
    ACTS_ERROR("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << state.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  // Post-stepping call to the action list
  ACTS_VERBOSE("Stepping loop done.");
  state.options.actionList(state, m_stepper, m_navigator, logger());

  // return progress flag here, decide on SUCCESS later
  return Result<void>::success();
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(const parameters_t& start,
                                       const propagator_options_t& options,
                                       bool makeCurvilinear) const
    -> Result<action_list_t_result_t<
        StepperCurvilinearTrackParameters,
        typename propagator_options_t::action_list_type>> {
  static_assert(
      std::is_copy_constructible<StepperCurvilinearTrackParameters>::value,
      "return track parameter type must be copy-constructible");

  auto state = makeState(start, options);

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, options,
                    makeCurvilinear);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(const parameters_t& start,
                                       const Surface& target,
                                       const propagator_options_t& options)
    const -> Result<action_list_t_result_t<
              StepperBoundTrackParameters,
              typename propagator_options_t::action_list_type>> {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  auto state = makeState<parameters_t, propagator_options_t, target_aborter_t,
                         path_aborter_t>(start, target, options);

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, target, options);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const parameters_t& start, const propagator_options_t& options) const {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperCurvilinearTrackParameters;

  static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                "return track parameter type must be copy-constructible");

  // Expand the abort list with a path aborter
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto abortList = options.abortList.append(pathAborter);

  // The expanded options (including path limit)
  auto eOptions = options.extend(abortList);
  eOptions.navigation.startSurface = &start.referenceSurface();
  eOptions.navigation.targetSurface = nullptr;
  using OptionsType = decltype(eOptions);
  using StateType =
      action_list_t_state_t<OptionsType,
                            typename propagator_options_t::action_list_type>;
  // Initialize the internal propagator state
  StateType state{
      eOptions,
      m_stepper.makeState(eOptions.geoContext, eOptions.magFieldContext, start,
                          eOptions.stepping.maxStepSize),
      m_navigator.makeState(eOptions.navigation)};

  static_assert(
      Concepts::has_method<const S, Result<double>, Concepts::Stepper::step_t,
                           StateType&, const N&>,
      "Step method of the Stepper is not compatible with the propagator "
      "state");

  initialize<StateType, path_aborter_t>(state);

  return state;
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  // Type of provided options
  target_aborter_t targetAborter;
  targetAborter.surface = &target;
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;
  auto abortList = options.abortList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(abortList);
  eOptions.navigation.startSurface = &start.referenceSurface();
  eOptions.navigation.targetSurface = &target;
  using OptionsType = decltype(eOptions);

  // Initialize the internal propagator state
  using StateType =
      action_list_t_state_t<OptionsType,
                            typename propagator_options_t::action_list_type>;
  StateType state{
      eOptions,
      m_stepper.makeState(eOptions.geoContext, eOptions.magFieldContext, start,
                          eOptions.stepping.maxStepSize),
      m_navigator.makeState(eOptions.navigation)};

  static_assert(
      Concepts::has_method<const S, Result<double>, Concepts::Stepper::step_t,
                           StateType&, const N&>,
      "Step method of the Stepper is not compatible with the propagator "
      "state");

  initialize<StateType, path_aborter_t>(state);

  return state;
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(propagator_state_t state,
                                        Result<void> propagationResult,
                                        const propagator_options_t& /*options*/,
                                        bool makeCurvilinear) const
    -> Result<action_list_t_result_t<
        StepperCurvilinearTrackParameters,
        typename propagator_options_t::action_list_type>> {
  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperCurvilinearTrackParameters;

  static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actions
  using ResultType =
      action_list_t_result_t<ReturnParameterType,
                             typename propagator_options_t::action_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  if (makeCurvilinear) {
    if (!m_stepper.prepareCurvilinearState(state, m_navigator)) {
      // information to compute curvilinearState is incomplete.
      return propagationResult.error();
    }
    /// Convert into return type and fill the result object
    auto curvState = m_stepper.curvilinearState(state.stepping);
    // Fill the end parameters
    result.endParameters =
        std::get<StepperCurvilinearTrackParameters>(curvState);
    // Only fill the transport jacobian when covariance transport was done
    if (state.stepping.covTransport) {
      result.transportJacobian = std::get<Jacobian>(curvState);
    }
  }

  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(propagator_state_t state,
                                        Result<void> propagationResult,
                                        const Surface& target,
                                        const propagator_options_t& /*options*/)
    const -> Result<action_list_t_result_t<
              StepperBoundTrackParameters,
              typename propagator_options_t::action_list_type>> {
  // Type of track parameters produced at the end of the propagation
  using ReturnParameterType = StepperBoundTrackParameters;

  static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actions
  using ResultType =
      action_list_t_result_t<ReturnParameterType,
                             typename propagator_options_t::action_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  // Compute the final results and mark the propagation as successful
  auto bsRes = m_stepper.boundState(state.stepping, target);
  if (!bsRes.ok()) {
    return bsRes.error();
  }
  const auto& bs = *bsRes;

  // Fill the end parameters
  result.endParameters = std::get<StepperBoundTrackParameters>(bs);
  // Only fill the transport jacobian when covariance transport was done
  if (state.stepping.covTransport) {
    result.transportJacobian = std::get<Jacobian>(bs);
  }
  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename path_aborter_t>
void Acts::Propagator<S, N>::initialize(propagator_state_t& state) const {
  state.position = m_stepper.position(state.stepping);
  state.direction = m_stepper.direction(state.stepping);

  // Navigator initialize state call
  m_navigator.initialize(state.navigation, state.position,
                         state.options.direction * state.direction);

  // Apply the loop protection - it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.abortList.template get<path_aborter_t>(),
      false, logger());
}

template <typename S, typename N>
template <typename propagator_state_t, typename result_t>
void Acts::Propagator<S, N>::moveStateToResult(propagator_state_t& state,
                                               result_t& result) const {
  result.tuple() = std::move(state.tuple());

  result.steps = state.steps;
  result.pathLength = state.pathLength;
}

template <typename derived_t>
Acts::Result<Acts::BoundTrackParameters>
Acts::detail::BasePropagatorHelper<derived_t>::propagateToSurface(
    const BoundTrackParameters& start, const Surface& target,
    const Options& options) const {
  using ResultType = Result<typename derived_t::template action_list_t_result_t<
      BoundTrackParameters, ActionList<>>>;
  using DerivedOptions = typename derived_t::template Options<>;

  DerivedOptions derivedOptions(options);

  // dummy initialization
  ResultType res = ResultType::failure(PropagatorError::Failure);

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

  if (res.ok()) {
    // Without errors we can expect a valid endParameters when propagating to a
    // target surface
    assert((*res).endParameters);
    return std::move((*res).endParameters.value());
  } else {
    return res.error();
  }
}
