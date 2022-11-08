// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// clang-format off
// Workaround for building on clang+libstdc++. Must be the first include.
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"
// clang-format on

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/VoidPropagatorComponents.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <functional>
#include <optional>
#include <type_traits>

#include <boost/algorithm/string.hpp>

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

  /// Final track parameters
  std::optional<parameters_t> endParameters = std::nullopt;

  /// Full transport jacobian
  std::optional<BoundMatrix> transportJacobian = std::nullopt;

  /// Number of propagation steps that were carried out
  unsigned int steps = 0;

  /// Signed distance over which the parameters were propagated
  double pathLength = 0.;
};

/// @brief Class holding the trivial options in propagator options
///
struct PropagatorPlainOptions {
  /// Propagation direction
  NavigationDirection direction = NavigationDirection::Forward;

  /// The |pdg| code for (eventual) material integration - pion default
  int absPdgCode = 211;

  /// The mass for the particle for (eventual) material integration
  double mass = 139.57018 * UnitConstants::MeV;

  /// Maximum number of steps for one propagate call
  unsigned int maxSteps = 1000;

  /// Maximum number of Runge-Kutta steps for the stepper step call
  unsigned int maxRungeKuttaStepTrials = 10000;

  /// Absolute maximum step size
  double maxStepSize = std::numeric_limits<double>::max();

  /// Absolute maximum path length
  double pathLimit = std::numeric_limits<double>::max();

  /// Required tolerance to reach target (surface, pathlength)
  double targetTolerance = s_onSurfaceTolerance;

  /// Loop protection step, it adapts the pathLimit
  bool loopProtection = true;
  double loopFraction = 0.5;  ///< Allowed loop fraction, 1 is a full loop

  // Configurations for Stepper
  /// Tolerance for the error of the integration
  double tolerance = 1e-4;

  /// Cut-off value for the step size
  double stepSizeCutOff = 0.;
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
struct PropagatorOptions : public PropagatorPlainOptions {
  using action_list_type = action_list_t;
  using aborter_list_type = aborter_list_t;

  /// Delete default contructor
  PropagatorOptions() = delete;

  /// PropagatorOptions copy constructor
  PropagatorOptions(
      const PropagatorOptions<action_list_t, aborter_list_t>& po) = default;

  /// PropagatorOptions with context
  PropagatorOptions(const GeometryContext& gctx,
                    const MagneticFieldContext& mctx, LoggerWrapper logger_)
      : geoContext(gctx), magFieldContext(mctx), logger(logger_) {}

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  PropagatorOptions<action_list_t, extended_aborter_list_t> extend(
      extended_aborter_list_t aborters) const {
    PropagatorOptions<action_list_t, extended_aborter_list_t> eoptions(
        geoContext, magFieldContext, logger);
    // Copy the options over
    eoptions.direction = direction;
    eoptions.absPdgCode = absPdgCode;
    eoptions.mass = mass;
    eoptions.maxSteps = maxSteps;
    eoptions.maxRungeKuttaStepTrials = maxRungeKuttaStepTrials;
    eoptions.maxStepSize = direction * std::abs(maxStepSize);
    eoptions.targetTolerance = targetTolerance;
    eoptions.pathLimit = direction * std::abs(pathLimit);
    eoptions.loopProtection = loopProtection;
    eoptions.loopFraction = loopFraction;

    // Stepper options
    eoptions.tolerance = tolerance;
    eoptions.stepSizeCutOff = stepSizeCutOff;
    // Action / abort list
    eoptions.actionList = std::move(actionList);
    eoptions.abortList = std::move(aborters);
    // And return the options
    return eoptions;
  }

  /// @brief Set the plain options
  ///
  /// @param pOptions The plain options
  void setPlainOptions(const PropagatorPlainOptions& pOptions) {
    // Copy the options over
    direction = pOptions.direction;
    absPdgCode = pOptions.absPdgCode;
    mass = pOptions.mass;
    maxSteps = pOptions.maxSteps;
    maxRungeKuttaStepTrials = pOptions.maxRungeKuttaStepTrials;
    maxStepSize = direction * std::abs(pOptions.maxStepSize);
    targetTolerance = pOptions.targetTolerance;
    pathLimit = direction * std::abs(pOptions.pathLimit);
    loopProtection = pOptions.loopProtection;
    loopFraction = pOptions.loopFraction;
    tolerance = pOptions.tolerance;
    stepSizeCutOff = pOptions.stepSizeCutOff;
  }

  /// List of actions
  action_list_t actionList;

  /// List of abort conditions
  aborter_list_t abortList;

  /// The context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// The context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  LoggerWrapper logger;
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
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;

  static_assert(StepperStateConcept<typename stepper_t::State>,
                "Stepper does not fulfill stepper concept.");
  static_assert(StepperConcept<stepper_t>,
                "Stepper does not fulfill stepper concept.");

 public:
  /// Type of the stepper in use for public scope
  using Stepper = stepper_t;

  /// Type of the navigator in use for public scope
  using Navigator = navigator_t;

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
    /// @param steppingIn Stepper state instance to begin with
    template <typename parameters_t>
    State(const parameters_t& start, const propagator_options_t& topts,
          StepperState steppingIn)
        : options(topts),
          stepping{std::move(steppingIn)},
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

 public:
  /// @brief Short-hand type definition for propagation result derived from
  ///        an action list
  ///
  /// @tparam parameters_t Type of the final track parameters
  /// @tparam action_list_t List of propagation action types
  ///
  template <typename parameters_t, typename action_list_t>
  using action_list_t_result_t =
      typename result_type_helper<parameters_t, action_list_t>::type;

 private:
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
  /// @tparam propagator_state_t Type of the propagator state with options
  ///
  /// @param [in,out] state the propagator state object
  /// @param [in,out] result an existing result object to start from
  ///
  /// @return Propagation result
  template <typename result_t, typename propagator_state_t>
  Result<void> propagate_impl(propagator_state_t& state,
                              result_t& result) const;

 public:
  /// @brief Propagate track parameters
  ///
  /// This function performs the propagation of the track parameters using the
  /// internal stepper implementation, until at least one abort condition is
  /// fulfilled or the maximum number of steps/path length provided in the
  /// propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] start initial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t, typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  Result<
      action_list_t_result_t<CurvilinearTrackParameters,
                             typename propagator_options_t::action_list_type>>
  propagate(const parameters_t& start,
            const propagator_options_t& options) const;

  /// @brief Propagate track parameters
  ///
  /// This function performs the propagation of the track parameters using the
  /// internal stepper implementation, until at least one abort condition is
  /// fulfilled or the maximum number of steps/path length provided in the
  /// propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] start initial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  /// @param [in] inputResult an existing result object to start from
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t, typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  Result<
      action_list_t_result_t<CurvilinearTrackParameters,
                             typename propagator_options_t::action_list_type>>
  propagate(
      const parameters_t& start, const propagator_options_t& options,
      action_list_t_result_t<CurvilinearTrackParameters,
                             typename propagator_options_t::action_list_type>&&
          inputResult) const;

  /// @brief Propagate track parameters - User method
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam target_aborter_t The target aborter type to be added
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  template <typename parameters_t, typename propagator_options_t,
            typename target_aborter_t = SurfaceReached,
            typename path_aborter_t = PathLimitReached>
  Result<action_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::action_list_type>>
  propagate(const parameters_t& start, const Surface& target,
            const propagator_options_t& options) const;

  /// @brief Propagate track parameters - User method
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @tparam parameters_t Type of initial track parameters to propagate
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam target_aborter_t The target aborter type to be added
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  /// @param [in] inputResult an existing result object to start from
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  template <typename parameters_t, typename propagator_options_t,
            typename target_aborter_t = SurfaceReached,
            typename path_aborter_t = PathLimitReached>
  Result<action_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::action_list_type>>
  propagate(
      const parameters_t& start, const Surface& target,
      const propagator_options_t& options,
      action_list_t_result_t<BoundTrackParameters,
                             typename propagator_options_t::action_list_type>
          inputResult) const;

 private:
  /// Implementation of propagation algorithm
  stepper_t m_stepper;

  /// Implementation of navigator
  navigator_t m_navigator;
};

}  // namespace Acts

#include "Acts/Propagator/Propagator.ipp"
