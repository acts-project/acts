// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <memory>
#include <type_traits>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/VoidPropagatorComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Result status of track parameter propagation
enum struct PropagatorStatus {
  SUCCESS,
  FAILURE,
  UNSET,
  IN_PROGRESS,
  WRONG_DIRECTION
};

/// @brief Simple class holding result of propagation call
///
/// @tparam parameters_t Type of final track parameters
/// @tparam result_list  Result pack for additional propagation
///                      quantities
template <typename parameters_t, typename... result_list>
struct Result : private detail::Extendable<result_list...>
{
  /// Constructor from initial propagation status
  ///
  /// @param s is the current status of the Result object
  Result(PropagatorStatus s = PropagatorStatus::UNSET)
    : detail::Extendable<result_list...>(), status(s)
  {
  }

  /// Accessor to additional propagation quantities
  using detail::Extendable<result_list...>::get;

  /// Final track parameters - initialized to null pointer
  std::unique_ptr<const parameters_t> endParameters = nullptr;

  /// Full transport jacobian
  std::unique_ptr<const ActsMatrixD<5, 5>> transportJacobian = nullptr;

  /// Propagation status
  PropagatorStatus status = PropagatorStatus::UNSET;

  /// Number of propagation steps that were carried out
  unsigned int steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;

  /// @brief Check the validity of the propagation result
  ///
  /// @return @c true if the final parameters are set and propagation status
  ///         is SUCCESS, otherwise @c false
  ///
  operator bool() const
  {
    return (endParameters && status == PropagatorStatus::SUCCESS);
  }
};

/// @brief Options for propagate() call
///
/// @tparam action_list_t List of action types called after each
///    propagation step with the current propagation and stepper state
///
/// @tparam aborter_list_t List of abort conditions tested after each
///    propagation step using the current propagation and stepper state
///
template <typename action_list_t  = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct PropagatorOptions
{

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  PropagatorOptions<action_list_t, extended_aborter_list_t>
  extend(extended_aborter_list_t aborters) const
  {
    PropagatorOptions<action_list_t, extended_aborter_list_t> eoptions;
    // Copy the options over
    eoptions.direction       = direction;
    eoptions.absPdgCode      = absPdgCode;
    eoptions.mass            = mass;
    eoptions.maxSteps        = maxSteps;
    eoptions.maxStepSize     = maxStepSize;
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit       = pathLimit;
    eoptions.loopProtection  = loopProtection;
    eoptions.loopFraction    = loopFraction;
    // Output option
    eoptions.debug         = debug;
    eoptions.debugString   = debugString;
    eoptions.debugPfxWidth = debugPfxWidth;
    eoptions.debugMsgWidth = debugMsgWidth;
    // Action / abort list
    eoptions.actionList = actionList;
    eoptions.abortList  = std::move(aborters);
    // And return the options
    return eoptions;
  }

  /// Propagation direction
  NavigationDirection direction = forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * units::_MeV;

  /// Maximum number of steps for one propagate() call
  unsigned int maxSteps = 1000;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Loop protection step, it adapts the pathLimit
  bool   loopProtection = true;
  double loopFraction   = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  /// Debug output steering:
  //  -> @todo: move to a debug struct
  // - the string where debug messages are stored (optionally)
  // - it also has some formatting options
  bool        debug         = false;  ///< switch debug on
  std::string debugString   = "";     ///< the string to collect msgs
  size_t      debugPfxWidth = 30;     ///< the prefix width
  size_t      debugMsgWidth = 50;     ///< the mesage width

  // Configurations for Stepper
  /// Tolerance for the error of the integration
  double tolerance = 1e-4;
  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;
};

/// @brief Propagator for particles (optionally in a magnetic field)
///
/// The Propagator works with a state objects given at function call
/// This state object contains the thread local state objects
///  - Navigator::state_type for object navigation and screen output
///  - Stepper::state_type state for the actual transport caching
///  (pos,dir,field)
///
/// @tparam stepper_t Type of stepper implementation of the propagation
/// @tparam naviagor_t Type of the navigator (optional)
///
/// This Propagator class serves as high-level steering code for propagating
/// track parameters. The actual implementation of the propagation has to be
/// implemented in the stepper_t object, which has to provide the following:
///
/// - a function for performing a single propagation step
/// - a type mapping for: initial track parameter type -> type of final track
///   parameters
/// - a type mapping for: (initial track parameter type and destination
///   surface type) -> type of final track parameters
/// - a type mapping for: initial track parameter type -> type of internal
///   state object
/// - a type mapping for: (initial track parameter type and destination
///   surface type) -> type of internal state object
///
template <typename stepper_t, typename navigator_t = detail::VoidNavigator>
class Propagator final
{
public:
  /// Type of the stepper in use for public scope
  using Stepper = stepper_t;

  /// Type of state object used by the propagation implementation
  using StepperState = typename stepper_t::template state_type<TrackParameters>;

  /// Typedef the navigator state
  using NavigatorState = typename navigator_t::state_type;

  /// Constructor from implementation object
  ///
  /// @param stepper The stepper implementation is moved to a private member
  /// @param navigator The navigator implementation, moved to a private member
  explicit Propagator(stepper_t stepper, navigator_t navigator = navigator_t())
    : m_stepper(std::move(stepper)), m_navigator(std::move(navigator))
  {
  }

  /// @brief private Propagator state for navigation and debugging
  ///
  /// @tparam parameters_t Type of the track parameters
  /// @tparam propagator_options_t Type of the Objections object
  ///
  /// This struct holds the common state information for propagating
  /// which is independent of the actual stepper implementation.
  template <typename propagator_options_t>
  struct State
  {

    /// Create the propagator state from the options
    ///
    /// @tparam parameters_t the type of the start parameters
    /// @tparam propagator_options_t the type of the propagator options
    ///
    /// @param start The start parameters, used to initialize stepping state
    /// @param topts The options handed over by the propagate call
    /// @param tabs The internal target aborters created in the call nethod
    template <typename parameters_t>
    State(const parameters_t& start, const propagator_options_t& topts)
      : options(topts), stepping(start, options.direction, options.maxStepSize)
    {
      // Setting the start surface
      navigation.startSurface = &start.referenceSurface();
    }

    /// These are the options - provided for each propagation step
    propagator_options_t options;

    /// Stepper state - internal state of the Stepper
    StepperState stepping;

    /// Navigation state - internal state of the Navigator
    NavigatorState navigation;
  };

private:
  /// @brief Helper struct determining the result's type
  ///
  /// @tparam parameters_t Type of final track parameters
  /// @tparam action_list_t    List of propagation action types
  ///
  /// This helper struct provides type definitions to extract the correct
  /// propagation result type from a given TrackParameter type and an
  /// ActionList.
  ///
  template <typename parameters_t, typename action_list_t>
  struct result_type_helper
  {
    /// @brief Propagation result type for an arbitrary list of additional
    ///        propagation results
    ///
    /// @tparam args Parameter pack specifying additional propagation results
    ///
    template <typename... args>
    using this_result_type = Result<parameters_t, args...>;

    /// @brief Propagation result type derived from a given action list
    using type = typename action_list_t::template result_type<this_result_type>;
  };

  /// @brief Short-hand type definition for propagation result derived from
  ///        an action list
  ///
  /// @tparam parameters_t Type of the final track parameters
  /// @tparam action_list_t List of propagation action types
  ///
  template <typename parameters_t, typename action_list_t>
  using action_list_t_result_t =
      typename result_type_helper<parameters_t, action_list_t>::type;

  /// @brief Propagate track parameters
  /// Private method with propagator and stepper state
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @note Does not (yet) convert into  the return_type of the propagation
  ///
  /// @tparam result_t Type of the result object for this propagation
  /// @tparam propagator_state_t Type of of propagator state with options
  ///
  /// @param [in,out] result of the propagation
  /// @param [in,out] state the propagator state object
  ///
  /// @return Propagation PropagatorStatus
  template <typename result_t, typename propagator_state_t>
  PropagatorStatus
  propagate_impl(result_t& result, propagator_state_t& state) const
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

public:
  /// @brief Propagate track parameters
  ///
  /// This function performs the propagation of the track parameters using the
  /// internal stepper implementation, until at least one abort condition is
  /// fulfilled or the maximum number of steps/path length provided in the
  /// propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam action_list_t Type list of actions, type ActionList<>
  /// @tparam aborter_list_t Type list of abort conditions, type AbortList<>
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start  nitial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t,
            typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename path_aborter_t = detail::PathLimitReached>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t>,
      action_list_t>
  propagate(
      const parameters_t& start,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced by the propagation
    using ReturnParameterType =
        typename stepper_t::template return_parameter_type<parameters_t>;

    // Type of the full propagation result, including output from actions
    using ResultType
        = action_list_t_result_t<ReturnParameterType, action_list_t>;

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

  /// @brief Propagate track parameters - User method
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam surface_t Type of target surface
  /// @tparam action_list_t Type list of actions
  /// @tparam aborter_list_t Type list of abort conditions
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  template <typename parameters_t,
            typename surface_t,
            typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename target_aborter_t = detail::SurfaceReached,
            typename path_aborter_t   = detail::PathLimitReached>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t,
                                                         surface_t>,
      action_list_t>
  propagate(
      const parameters_t& start,
      const surface_t&    target,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced at the end of the propagation
    using return_parameter_type =
        typename stepper_t::template return_parameter_type<parameters_t,
                                                           surface_t>;

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

private:
  /// Implementation of propagation algorithm
  stepper_t m_stepper;

  /// Implementation of navigator
  navigator_t m_navigator;

  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// options.debug == true case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the nested propagator state object
  ///
  /// @param state the propagator state for the debug flag, prefix/length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
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
};

}  // namespace Acts
