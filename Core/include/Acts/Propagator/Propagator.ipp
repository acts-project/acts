// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"

#include <type_traits>

template <typename S, typename N>
template <typename result_t, typename propagator_state_t>
auto Acts::Propagator<S, N>::propagate_impl(propagator_state_t& state,
                                            result_t& result) const
    -> Result<void> {
  // Pre-stepping call to the navigator and action list
  ACTS_VERBOSE("Entering propagation.");

  // Navigator initialize state call
  m_navigator.initialize(state, m_stepper);
  // Pre-Stepping call to the action list
  state.options.actionList(state, m_stepper, m_navigator, result, logger());
  // assume negative outcome, only set to true later if we actually have
  // a positive outcome.

  // start at true, if we don't begin the stepping loop we're fine.
  bool terminatedNormally = true;

  // Pre-Stepping: abort condition check
  if (!state.options.abortList(state, m_stepper, m_navigator, result,
                               logger())) {
    // Stepping loop
    ACTS_VERBOSE("Starting stepping loop.");

    terminatedNormally = false;  // priming error condition

    // Propagation loop : stepping
    for (; result.steps < state.options.maxSteps; ++result.steps) {
      // Pre-Stepping: target setting
      m_navigator.preStep(state, m_stepper);
      // Perform a propagation step - it takes the propagation state
      Result<double> res = m_stepper.step(state, m_navigator);
      if (res.ok()) {
        // Accumulate the path length
        double s = *res;
        result.pathLength += s;
        ACTS_VERBOSE("Step with size = " << s << " performed");
      } else {
        ACTS_ERROR("Step failed with " << res.error() << ": "
                                       << res.error().message());
        // pass error to caller
        return res.error();
      }
      // Post-stepping:
      // navigator post step call - action list - aborter list
      m_navigator.postStep(state, m_stepper);
      state.options.actionList(state, m_stepper, m_navigator, result, logger());
      if (state.options.abortList(state, m_stepper, m_navigator, result,
                                  logger())) {
        terminatedNormally = true;
        break;
      }
    }
  } else {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");
  }

  // if we didn't terminate normally (via aborters) set navigation break.
  // this will trigger error output in the lines below
  if (!terminatedNormally) {
    m_navigator.navigationBreak(state.navigation, true);
    ACTS_ERROR("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << result.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  // Post-stepping call to the action list
  ACTS_VERBOSE("Stepping loop done.");
  state.options.actionList(state, m_stepper, m_navigator, result, logger());

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
  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperCurvilinearTrackParameters;

  static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actions
  using ResultType =
      action_list_t_result_t<ReturnParameterType,
                             typename propagator_options_t::action_list_type>;

  return propagate<parameters_t, propagator_options_t, path_aborter_t>(
      start, options, makeCurvilinear, ResultType{});
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(
    const parameters_t& start, const propagator_options_t& options,
    bool makeCurvilinear,
    action_list_t_result_t<StepperCurvilinearTrackParameters,
                           typename propagator_options_t::action_list_type>&&
        inputResult) const
    -> Result<action_list_t_result_t<
        StepperCurvilinearTrackParameters,
        typename propagator_options_t::action_list_type>> {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  using ResultType = std::decay_t<decltype(inputResult)>;

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
  using OptionsType = decltype(eOptions);
  // Initialize the internal propagator state
  using StateType = State<OptionsType>;
  StateType state{
      eOptions,
      m_stepper.makeState(eOptions.geoContext, eOptions.magFieldContext, start,
                          eOptions.maxStepSize),
      m_navigator.makeState(&start.referenceSurface(), nullptr)};

  static_assert(
      Concepts::has_method<const S, Result<double>, Concepts::Stepper::step_t,
                           StateType&, const N&>,
      "Step method of the Stepper is not compatible with the propagator "
      "state");

  // Apply the loop protection - it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.abortList.template get<path_aborter_t>(),
      false, logger());
  // Perform the actual propagation & check its outcome
  auto result = propagate_impl<ResultType>(state, inputResult);
  if (result.ok()) {
    if (makeCurvilinear) {
      /// Convert into return type and fill the result object
      auto curvState = m_stepper.curvilinearState(state.stepping);
      // Fill the end parameters
      inputResult.endParameters =
          std::get<StepperCurvilinearTrackParameters>(curvState);
      // Only fill the transport jacobian when covariance transport was done
      if (state.stepping.covTransport) {
        inputResult.transportJacobian = std::get<Jacobian>(curvState);
      }
    }
    return Result<ResultType>::success(std::forward<ResultType>(inputResult));
  } else {
    return result.error();
  }
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<action_list_t_result_t<
        StepperBoundTrackParameters,
        typename propagator_options_t::action_list_type>> {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  // Type of track parameters produced at the end of the propagation
  using return_parameter_type = StepperBoundTrackParameters;

  // Type of the full propagation result, including output from actions
  using ResultType =
      action_list_t_result_t<return_parameter_type,
                             typename propagator_options_t::action_list_type>;

  return propagate<parameters_t, propagator_options_t, target_aborter_t,
                   path_aborter_t>(start, target, options, ResultType{});
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options,
    action_list_t_result_t<StepperBoundTrackParameters,
                           typename propagator_options_t::action_list_type>
        inputResult) const
    -> Result<action_list_t_result_t<
        StepperBoundTrackParameters,
        typename propagator_options_t::action_list_type>> {
  static_assert(Concepts::BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  using ResultType = std::decay_t<decltype(inputResult)>;

  // Type of provided options
  target_aborter_t targetAborter;
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;
  auto abortList = options.abortList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(abortList);
  using OptionsType = decltype(eOptions);

  // Initialize the internal propagator state
  using StateType = State<OptionsType>;
  StateType state{
      eOptions,
      m_stepper.makeState(eOptions.geoContext, eOptions.magFieldContext, start,
                          eOptions.maxStepSize),
      m_navigator.makeState(&start.referenceSurface(), &target)};

  static_assert(
      Concepts::has_method<const S, Result<double>, Concepts::Stepper::step_t,
                           StateType&, const N&>,
      "Step method of the Stepper is not compatible with the propagator "
      "state");

  // Apply the loop protection, it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.abortList.template get<path_aborter_t>(),
      false, logger());

  // Perform the actual propagation
  auto result = propagate_impl<ResultType>(state, inputResult);

  if (result.ok()) {
    // Compute the final results and mark the propagation as successful
    auto bsRes = m_stepper.boundState(state.stepping, target);
    if (!bsRes.ok()) {
      return bsRes.error();
    }

    const auto& bs = *bsRes;

    // Fill the end parameters
    inputResult.endParameters = std::get<StepperBoundTrackParameters>(bs);
    // Only fill the transport jacobian when covariance transport was done
    if (state.stepping.covTransport) {
      inputResult.transportJacobian = std::get<Jacobian>(bs);
    }
    return Result<ResultType>::success(std::forward<ResultType>(inputResult));
  } else {
    return result.error();
  }
}
