// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename S, typename N>
template <typename result_t, typename propagator_state_t>
Acts::PropagatorStatus
Acts::Propagator<S, N>::propagate_impl(result_t& result,
                                       propagator_state_t& state) const
{

  // Pre-stepping call to the navigator and action list
  debugLog(state, [&] { return std::string("Entering propagation."); });

  // Navigator initialize state call
  m_navigator.status(state, m_stepper);
  // Pre-Stepping call to the action list
  state.options.actionList(state, m_stepper, result);
  // assume negative outcome, only set to true later if we actually have
  // a positive outcome.
  // This is needed for correct error logging
  bool terminatedNormally = false;
  // Pre-Stepping: abort condition check
  if (!state.options.abortList(result, state, m_stepper)) {
    // Pre-Stepping: target setting
    m_navigator.target(state, m_stepper);
    // Stepping loop
    debugLog(state, [&] { return std::string("Starting stepping loop."); });
    // Propagation loop : stepping
    for (; result.steps < state.options.maxSteps; ++result.steps) {
      // Perform a propagation step - it takes the propagation state
      double s = m_stepper.step(state);
      // Accumulate the path length
      result.pathLength += s;
      // Call the actions, can (& will likely) modify the state
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Step with size = ";
        dstream << s;
        dstream << " performed.";
        return dstream.str();
      });
      // Post-step
      // navigator status call - action list - aborter list - target call
      m_navigator.status(state, m_stepper);
      state.options.actionList(state, m_stepper, result);
      if (state.options.abortList(result, state, m_stepper)) {
        terminatedNormally = true;
        break;
      }
      m_navigator.target(state, m_stepper);
    }
  }

  // if we didn't terminate normally (via aborters) set navigation break.
  // this will trigger error output in the lines below
  if (!terminatedNormally) {
    state.navigation.navigationBreak = true;
  }

  // Post-stepping call to the action list
  debugLog(state, [&] { return std::string("Stepping loop done."); });
  state.options.actionList(state, m_stepper, result);

  // return progress flag here, decide on SUCCESS later
  return PropagatorStatus::IN_PROGRESS;
}

template <typename S, typename N>
template <typename parameters_t,
          typename action_list_t,
          typename aborter_list_t,
          template <typename, typename> class propagator_options_t,
          typename path_aborter_t>
auto
Acts::Propagator<S, N>::propagate(
    const parameters_t& start,
    const propagator_options_t<action_list_t, aborter_list_t>& options) const
    -> action_list_t_result_t<
        typename S::template return_parameter_type<parameters_t>,
        action_list_t>
{

  // Type of track parameters produced by the propagation
  using ReturnParameterType =
      typename S::template return_parameter_type<parameters_t>;

  // Type of the full propagation result, including output from actions
  using ResultType = action_list_t_result_t<ReturnParameterType, action_list_t>;

  static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                "return track parameter type must be copy-constructible");

  // Initialize the propagation result object
  ResultType result(PropagatorStatus::IN_PROGRESS);

  // Expand the abort list with a path aborter
  path_aborter_t pathAborter;
  auto           abortList = options.abortList.append(pathAborter);

  // The expanded options (including path limit)
  auto eOptions     = options.extend(abortList);
  using OptionsType = decltype(eOptions);
  // Initialize the internal propagator state
  using StateType = State<OptionsType>;
  StateType state(start, eOptions);

  // Apply the loop protection - it resets the internal path limit
  if (options.loopProtection) {
    detail::LoopProtection<path_aborter_t> lProtection;
    lProtection(state, m_stepper);
  }

  // Perform the actual propagation & check its outcome
  if (propagate_impl(result, state) != PropagatorStatus::IN_PROGRESS) {
    result.status = PropagatorStatus::FAILURE;
  } else {
    /// Convert into return type and fill the result object
    m_stepper.convert(state.stepping, result);
    result.status = PropagatorStatus::SUCCESS;
  }

  return result;
}

template <typename S, typename N>
template <typename parameters_t,
          typename surface_t,
          typename action_list_t,
          typename aborter_list_t,
          template <typename, typename> class propagator_options_t,
          typename target_aborter_t,
          typename path_aborter_t>
auto
Acts::Propagator<S, N>::propagate(
    const parameters_t& start,
    const surface_t&    target,
    const propagator_options_t<action_list_t, aborter_list_t>& options) const
    -> action_list_t_result_t<
        typename S::template return_parameter_type<parameters_t, surface_t>,
        action_list_t>
{

  // Type of track parameters produced at the end of the propagation
  using return_parameter_type =
      typename S::template return_parameter_type<parameters_t, surface_t>;

  // Type of provided options
  target_aborter_t targetAborter;
  path_aborter_t   pathAborter;
  auto abortList = options.abortList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions     = options.extend(abortList);
  using OptionsType = decltype(eOptions);

  // Type of the full propagation result, including output from actions
  using ResultType
      = action_list_t_result_t<return_parameter_type, action_list_t>;

  // Initialize the propagation result object
  ResultType result(PropagatorStatus::IN_PROGRESS);

  static_assert(std::is_copy_constructible<return_parameter_type>::value,
                "return track parameter type must be copy-constructible");

  // Initialize the internal propagator state
  using StateType = State<OptionsType>;
  StateType state(start, eOptions);
  state.navigation.targetSurface = &target;

  // Apply the loop protection, it resets the interal path limit
  detail::LoopProtection<path_aborter_t> lProtection;
  lProtection(state, m_stepper);

  // Perform the actual propagation
  if (propagate_impl(result, state) != PropagatorStatus::IN_PROGRESS) {
    result.status = PropagatorStatus::FAILURE;
  } else {
    // Compute the final results and mark the propagation as successful
    m_stepper.convert(state.stepping, result, target);
    result.status = PropagatorStatus::SUCCESS;
  }
  return result;
}

template <typename S, typename N>
template <typename propagator_state_t>
void
Acts::Propagator<S, N>::debugLog(
    propagator_state_t&                 state,
    const std::function<std::string()>& logAction) const
{
  if (state.options.debug) {
    std::vector<std::string> lines;
    std::string              input = logAction();
    boost::split(lines, input, boost::is_any_of("\n"));
    for (const auto& line : lines) {
      std::stringstream dstream;
      dstream << "|->" << std::setw(state.options.debugPfxWidth);
      dstream << "Propagator"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
      state.options.debugString += dstream.str();
    }
  }
}
