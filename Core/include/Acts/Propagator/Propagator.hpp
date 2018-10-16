// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <memory>
#include <type_traits>
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAbortConditions.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Result status of track parameter propagation
enum struct Status { SUCCESS, FAILURE, UNSET, IN_PROGRESS, WRONG_DIRECTION };

/// @brief The void navigator struct as a default navigator
///
/// It does not provide any navigation action, the compiler
/// should eventually optimise that the function core is not done
struct VoidNavigator
{

  /// Nested State struct
  struct State
  {
    /// Navigation state - external state: the start surface
    const Surface* startSurface = nullptr;

    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation state - external state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator if the target is reached
    bool targetReached = false;

    /// Navigation state : a break has been detected
    bool navigationBreak = false;
  };

  /// Unique typedef to publish to the Propagator
  using state_type = State;

  /// Navigation call
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  ///
  /// Empty call, hopefully the compiler checks this
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& /*state*/) const
  {
  }
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
  Result(Status s = Status::UNSET)
    : detail::Extendable<result_list...>(), status(s)
  {
  }

  /// Accessor to additional propagation quantities
  using detail::Extendable<result_list...>::get;

  /// Final track parameters - initialized to null pointer
  std::unique_ptr<const parameters_t> endParameters = nullptr;

  /// Propagation status
  Status status = Status::UNSET;

  /// Number of propagation steps that were carried out
  unsigned int steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;

  /// @brief Check the validity of the propagation result
  ///
  /// @return @c true if the final parameters are set and propagation status
  ///         is SUCCESS, otherwise @c false
  ///
  operator bool() const { return (endParameters && status == Status::SUCCESS); }
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

  /// Propagation direction
  NavigationDirection direction = forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * units::_MeV;

  /// Maximum number of steps for one propagate() call
  unsigned int maxSteps = 1000;

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Absolute maximum step size
  double maxStepSize = 1 * units::_m;

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Loop protection step, it adapts the pathLimit
  bool   loopProtection = true;
  double loopFraction   = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  /// Debug output steering:
  // - the string where debug messages are stored (optionally)
  // - it also has some formatting options
  bool        debug         = false;  ///< switch debug on
  std::string debugString   = "";     ///< the string to collect msgs
  size_t      debugPfxWidth = 30;     ///< the prefix width
  size_t      debugMsgWidth = 50;     ///< the mesage width

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t stopConditions;
};

/// @brief Propagator for particles (optionally in a magnetic field)
///
/// The Propagator works with a state objects given at function call
/// This state object contains the thread local state objects
///  - Navigator::state_type for object navigation and screen output
///  - Stepper::state_type state for the actual transport caching
///  (pos,dir,field)
///
/// @tparam stepper_t stepper implementation of the propagation algorithm
/// @tparam navigator_list_t the (optional) navigator type, it is a type
///         of action_list with is called before all the other optios
/// @tparam path_aborter_t Type of object to abort when path limit reached
/// @tparam target_aborter_t Type of the object to do initial step estimation
///                          and target abort
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
template <typename stepper_t, typename navigator_t = VoidNavigator>
class Propagator final
{
public:
  /// Type of state object used by the propagation implementation
  using StepperState = typename stepper_t::template state_type<TrackParameters>;

  /// Typedef the navigator state
  using NavigatorState = typename navigator_t::state_type;

  /// Constructor from implementation object
  ///
  /// @param stepper The stepper implementation is moved to a private member
  /// @param navigator The navigator implentiation, moved to a private member
  explicit Propagator(stepper_t stepper, navigator_t navigator = navigator_t())
    : m_stepper(std::move(stepper)), m_navigator(std::move(navigator))
  {
  }

private:
  /// @brief private Propagator state for navigation and debugging
  ///
  /// @tparam parameters_t Type of the track parameters
  /// @tparam propagator_options_t Type of the Objections object
  /// @tparam target_aborter_list_t Type of the aborter list
  ///
  /// This struct holds the common state information for propagating
  /// which is independent of the actual stepper implementation.
  template <typename parameters_t,
            typename propagator_option_t,
            typename target_aborter_list_t>
  struct State
  {

    /// Create the propagator state from the options
    ///
    /// @tparam parameters_t the type of the start parameters
    /// @tparam propagator_option_t the type of the propagator options
    /// @tparam target_aborter_list_t the type of the target aborters
    ///
    /// @param start The start parameters, used to initialize stepping state
    /// @param topts The options handed over by the propagate call
    /// @param tabs The internal target aborters created in the call nethod
    State(const parameters_t&        start,
          const propagator_option_t& topts,
          target_aborter_list_t      tabs)
      : options(topts)
      , targetAborters(std::move(tabs))
      , stepping(start, options.direction, options.maxStepSize)
    {
    }

    /// These are the options - provided for each propagation step
    propagator_option_t options;

    /// These are the target aborters (internally created)
    target_aborter_list_t targetAborters;

    /// Stepper state - internal state of the Stepper
    StepperState stepping;

    /// Navigation state - internal state of the Navigator
    NavigatorState navigation;
  };

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
  /// @tparam T       Type of the final track parameters
  /// @tparam action_list_t List of propagation action types
  ///
  template <typename T, typename action_list_t>
  using action_list_t_result_t =
      typename result_type_helper<T, action_list_t>::type;

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
  /// @return Propagation Status
  template <typename result_t, typename propagator_state_t>
  Status
  propagate_(result_t& result, propagator_state_t& state) const
  {
    // Pre-stepping call to the abort list
    debugLog(state,
             [&] { return std::string("Calling pre-stepping aborters."); });
    if (state.targetAborters(result, state)) {
      return Status::FAILURE;
    }
    // Pre-stepping call to the navigator and action list
    debugLog(state, [&] {
      return std::string("Calling pre-stepping navigator & actions.");
    });
    m_navigator(state);
    state.options.actionList(state, result);

    // Propagation loop : stepping
    for (; result.steps < state.options.maxSteps; ++result.steps) {
      // Perform a propagation step - it only takes the stepping state
      double s = m_stepper.step(state.stepping);
      // accumulate the path length
      result.pathLength += s;
      // Call the actions, can (& will likely) modify the state
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Calling navigator & actions after step of size = ";
        dstream << s;
        return dstream.str();
      });
      m_navigator(state);
      state.options.actionList(state, result);
      // Call the stop_conditions and the internal stop conditions
      // break condition triggered, but still count the step
      debugLog(state,
               [&] { return std::string("Calling aborters after step."); });
      if (state.options.stopConditions(result, state)
          || state.targetAborters(result, state)) {
        break;
      }
    }
    // Post-stepping call to the action list
    debugLog(state,
             [&] { return std::string("Calling post-stepping action list."); });

    state.options.actionList(state, result);

    // return progress flag here, decide on SUCCESS later
    return Status::IN_PROGRESS;
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
            typename path_arborter_t = detail::PathLimitReached>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t>,
      action_list_t>
  propagate(
      const parameters_t& start,
      const PropagatorOptions<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced by the propagation
    using ReturnParameterType =
        typename stepper_t::template return_parameter_type<parameters_t>;

    // Type of the full propagation result, including output from actions
    using ResultType
        = action_list_t_result_t<ReturnParameterType, action_list_t>;

    // Type of provided options which consist action and abort list
    using OptionsType = PropagatorOptions<action_list_t, aborter_list_t>;

    static_assert(std::is_copy_constructible<ReturnParameterType>::value,
                  "return track parameter type must be copy-constructible");

    // Initialize the propagation result object
    ResultType result(Status::IN_PROGRESS);

    // Internal Abort list - only with path limit as no target surface given
    using TargetAborters = AbortList<path_arborter_t>;
    TargetAborters targetAborters;

    // Initialize the internal propagator state
    using StateType = State<parameters_t, OptionsType, TargetAborters>;
    StateType state(start, options, targetAborters);

    // Apply the loop protection - it resets the internal path limit
    if (options.loopProtection) {
      detail::LoopProtection<path_arborter_t> lProtection;
      lProtection(state, m_stepper);
    }

    // Perform the actual propagation & check its outcome
    if (propagate_(result, state) != Status::IN_PROGRESS) {
      result.status = Status::FAILURE;
    } else {
      /// Convert into the return type
      result.endParameters = std::make_unique<const ReturnParameterType>(
          m_stepper.convert(state.stepping));
      result.status = Status::SUCCESS;
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
            typename path_arborter_t  = detail::PathLimitReached,
            typename target_aborter_t = detail::SurfaceReached>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t,
                                                         surface_t>,
      action_list_t>
  propagate(
      const parameters_t& start,
      const surface_t&    target,
      const PropagatorOptions<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced at the end of the propagation
    using return_parameter_type =
        typename stepper_t::template return_parameter_type<parameters_t,
                                                           surface_t>;

    // Type of provided options
    using OptionsType = PropagatorOptions<action_list_t, aborter_list_t>;

    // Type of the full propagation result, including output from actions
    using ResultType
        = action_list_t_result_t<return_parameter_type, action_list_t>;

    // Initialize the propagation result object
    ResultType result(Status::IN_PROGRESS);

    static_assert(std::is_copy_constructible<return_parameter_type>::value,
                  "return track parameter type must be copy-constructible");

    // Internal Abort list for target and path surface
    using TargetAborters = AbortList<target_aborter_t, path_arborter_t>;
    TargetAborters targetAborters;

    // Initialize the internal propagator state
    using StateType = State<parameters_t, OptionsType, TargetAborters>;
    StateType state(start, options, targetAborters);

    // Setting the start and the target surface
    state.navigation.startSurface  = &start.referenceSurface();
    state.navigation.targetSurface = &target;

    // Apply the loop protection, it resets the interal path limit
    detail::LoopProtection<path_arborter_t> lProtection;
    lProtection(state, m_stepper);

    // Perform the actual propagation
    if (propagate_(result, state) != Status::IN_PROGRESS) {
      result.status = Status::FAILURE;
    } else {
      // Compute the final results and mark the propagation as successful
      result.endParameters = std::make_unique<const return_parameter_type>(
          m_stepper.convert(state.stepping, target));
      result.status = Status::SUCCESS;
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
      std::stringstream dstream;
      dstream << "|->" << std::setw(state.options.debugPfxWidth);
      dstream << "Propagator"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};

}  // namespace Acts
