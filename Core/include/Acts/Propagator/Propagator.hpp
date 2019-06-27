// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <type_traits>

#include <boost/algorithm/string.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/VoidPropagatorComponents.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @brief Simple class holding result of propagation call
///
/// @tparam parameters_t Type of final track parameters
/// @tparam result_list  Result pack for additional propagation
///                      quantities
template <typename parameters_t, typename... result_list>
struct PropagatorResult : private detail::Extendable<result_list...> {
  /// Accessor to additional propagation quantities
  using detail::Extendable<result_list...>::get;

  /// Final track parameters - initialized to null pointer
  std::unique_ptr<const parameters_t> endParameters = nullptr;

  /// Full transport jacobian
  std::unique_ptr<const BoundMatrix> transportJacobian = nullptr;

  /// Number of propagation steps that were carried out
  unsigned int steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;
};

/// @brief Options for propagate() call
///
/// @tparam action_list_t List of action types called after each
///    propagation step with the current propagation and stepper state
///
/// @tparam aborter_list_t List of abort conditions tested after each
///    propagation step using the current propagation and stepper state
///
template <typename action_list_t = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct PropagatorOptions {
  /// Delete default contructor
  PropagatorOptions() = delete;

  /// PropagatorOptions copy constructor
  PropagatorOptions(
      const PropagatorOptions<action_list_t, aborter_list_t>& po) = default;

  /// PropagatorOptions with context
  PropagatorOptions(std::reference_wrapper<const GeometryContext> gctx,
                    std::reference_wrapper<const MagneticFieldContext> mctx)
      : geoContext(gctx), magFieldContext(mctx) {}

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  PropagatorOptions<action_list_t, extended_aborter_list_t> extend(
      extended_aborter_list_t aborters) const {
    PropagatorOptions<action_list_t, extended_aborter_list_t> eoptions(
        geoContext, magFieldContext);
    // Copy the options over
    eoptions.direction = direction;
    eoptions.absPdgCode = absPdgCode;
    eoptions.mass = mass;
    eoptions.maxSteps = maxSteps;
    eoptions.maxStepSize = maxStepSize;
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit = pathLimit;
    eoptions.loopProtection = loopProtection;
    eoptions.loopFraction = loopFraction;
    // Output option
    eoptions.debug = debug;
    eoptions.debugString = debugString;
    eoptions.debugPfxWidth = debugPfxWidth;
    eoptions.debugMsgWidth = debugMsgWidth;
    // Stepper options
    eoptions.tolerance = tolerance;
    eoptions.stepSizeCutOff = stepSizeCutOff;
    // Action / abort list
    eoptions.actionList = std::move(actionList);
    eoptions.abortList = std::move(aborters);
    // And return the options
    return eoptions;
  }

  /// Propagation direction
  NavigationDirection direction = forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * UnitConstants::MeV;

  /// Maximum number of steps for one propagate call
  unsigned int maxSteps = 1000;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Loop protection step, it adapts the pathLimit
  bool loopProtection = true;
  double loopFraction = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  /// Debug output steering:
  //  -> @todo: move to a debug struct
  // - the string where debug messages are stored (optionally)
  // - it also has some formatting options
  bool debug = false;            ///< switch debug on
  std::string debugString = "";  ///< the string to collect msgs
  size_t debugPfxWidth = 30;     ///< the prefix width
  size_t debugMsgWidth = 50;     ///< the mesage width

  // Configurations for Stepper
  /// Tolerance for the error of the integration
  double tolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
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
class Propagator final {
  using Jacobian = BoundMatrix;
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

  static_assert(StepperStateConcept<typename stepper_t::State>,
                "Stepper does not fulfill stepper concept.");
  static_assert(StepperConcept<stepper_t>,
                "Stepper does not fulfill stepper concept.");

 public:
  /// Type of the stepper in use for public scope
  using Stepper = stepper_t;

  /// Type of state object used by the propagation implementation
  using StepperState = typename Stepper::State;

  /// Typedef the navigator state
  using NavigatorState = typename navigator_t::State;

  /// Constructor from implementation object
  ///
  /// @param stepper The stepper implementation is moved to a private member
  /// @param navigator The navigator implementation, moved to a private member
  explicit Propagator(stepper_t stepper, navigator_t navigator = navigator_t())
      : m_stepper(std::move(stepper)), m_navigator(std::move(navigator)) {}

  /// @brief private Propagator state for navigation and debugging
  ///
  /// @tparam parameters_t Type of the track parameters
  /// @tparam propagator_options_t Type of the Objections object
  ///
  /// This struct holds the common state information for propagating
  /// which is independent of the actual stepper implementation.
  template <typename propagator_options_t>
  struct State {
    /// Create the propagator state from the options
    ///
    /// @tparam parameters_t the type of the start parameters
    /// @tparam propagator_options_t the type of the propagator options
    ///
    /// @param start The start parameters, used to initialize stepping state
    /// @param topts The options handed over by the propagate call
    template <typename parameters_t>
    State(const parameters_t& start, const propagator_options_t& topts)
        : options(topts),
          stepping(topts.geoContext, topts.magFieldContext, start,
                   topts.direction, topts.maxStepSize),
          geoContext(topts.geoContext) {
      // Setting the start surface
      navigation.startSurface = &start.referenceSurface();
    }

    /// These are the options - provided for each propagation step
    propagator_options_t options;

    /// Stepper state - internal state of the Stepper
    StepperState stepping;

    /// Navigation state - internal state of the Navigator
    NavigatorState navigation;

    /// Context object for the geometry
    std::reference_wrapper<const GeometryContext> geoContext;
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
  struct result_type_helper {
    /// @brief Propagation result type for an arbitrary list of additional
    ///        propagation results
    ///
    /// @tparam args Parameter pack specifying additional propagation results
    ///
    template <typename... args>
    using this_result_type = PropagatorResult<parameters_t, args...>;

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
  Result<result_t> propagate_impl(propagator_state_t& state) const;

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
  /// @param [in] start initial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t, typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename path_aborter_t = detail::PathLimitReached>
  Result<action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t>,
      action_list_t>>
  propagate(
      const parameters_t& start,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const;

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
  template <typename parameters_t, typename surface_t, typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t,
            typename target_aborter_t = detail::SurfaceReached,
            typename path_aborter_t = detail::PathLimitReached>
  Result<
      action_list_t_result_t<typename stepper_t::template return_parameter_type<
                                 parameters_t, surface_t>,
                             action_list_t>>
  propagate(
      const parameters_t& start, const surface_t& target,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const;

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
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const;
};

}  // namespace Acts

#include "Acts/Propagator/Propagator.ipp"