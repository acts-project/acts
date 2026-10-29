// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DefaultActor.hpp"
#include "Acts/Propagator/NavigatorConcept.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/PropagatorResult.hpp"
#include "Acts/Propagator/PropagatorState.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

/// Common simplified base interface for propagators.
///
/// This class only supports propagation from start bound parameters to a target
/// surface and returns only the end bound parameters.
/// Navigation is performed if the underlying propagator is configured with an
/// appropriate navigator. No custom actors or aborters are supported.
class BasePropagator {
 public:
  /// Base propagator options
  using Options = PropagatorPlainOptions;

  /// Method to propagate start bound track parameters to a target surface.
  /// @param start The start bound track parameters.
  /// @param target The target surface.
  /// @param options The propagation options.
  /// @return The end bound track parameters.
  virtual Result<BoundTrackParameters> propagateToSurface(
      const BoundTrackParameters& start, const Surface& target,
      const Options& options) const = 0;

  virtual ~BasePropagator() = default;
};

namespace detail {
class PropagatorStub {};

template <typename derived_t>
class BasePropagatorHelper : public BasePropagator {
 public:
  Result<BoundTrackParameters> propagateToSurface(
      const BoundTrackParameters& start, const Surface& target,
      const Options& options) const override;
};
}  // namespace detail

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
template <StepperConcept stepper_t,
          NavigatorConcept navigator_t = VoidNavigator>
class Propagator final
    : public std::conditional_t<
          SupportsBoundParameters_v<stepper_t>,
          detail::BasePropagatorHelper<Propagator<stepper_t, navigator_t>>,
          detail::PropagatorStub> {
  using BoundParameters = stepper_t::BoundParameters;

  using Jacobian = BoundMatrix;
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;

 public:
  /// Type of the stepper in use
  using Stepper = stepper_t;

  /// Type of the navigator in use
  using Navigator = navigator_t;

  /// Type of the stepper options
  using StepperOptions = typename Stepper::Options;

  /// Type of the navigator options
  using NavigatorOptions = typename Navigator::Options;

  /// Type of the propagator options with an actor
  template <typename actor_t = DefaultActor>
  using Options = PropagatorOptions<StepperOptions, NavigatorOptions, actor_t>;

  /// Type of the stepper state
  using StepperState = typename Stepper::State;

  /// Type of the navigator state
  using NavigatorState = typename Navigator::State;

  /// Type of the propagation state derived from the propagator options
  /// @tparam propagator_options_t Type of the propagator options
  template <typename propagator_options_t>
  using State =
      PropagatorState<propagator_options_t, StepperState, NavigatorState>;

  /// Type of the propagation result
  using ResultType = PropagatorResult<BoundParameters>;

  /// Constructor from implementation object
  ///
  /// @param stepper The stepper implementation is moved to a private member
  /// @param navigator The navigator implementation, moved to a private member
  /// @param _logger a logger instance
  explicit Propagator(Stepper stepper, Navigator navigator = Navigator(),
                      std::shared_ptr<const Logger> _logger =
                          getDefaultLogger("Propagator", Acts::Logging::INFO))
      : m_stepper(std::move(stepper)),
        m_navigator(std::move(navigator)),
        m_logger{std::move(_logger)} {}

  /// Access to the stepper instance
  /// @return Const reference to the stepper
  const Stepper& stepper() const { return m_stepper; }

  /// Access to the navigator instance
  /// @return Const reference to the navigator
  const Navigator& navigator() const { return m_navigator; }

  /// @brief Propagate track parameters
  ///
  /// This function propagates the track parameters until the actor of the
  /// options aborts or the maximum number of steps is reached.
  ///
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start initial track parameters to propagate
  /// @param [in] options Propagation options, type Options<,>
  /// @param [in] createFinalParameters Whether to produce parameters at the end of the propagation
  ///
  /// @return Propagation result containing the propagation status and the
  ///         final track parameters
  template <typename propagator_options_t>
  Result<ResultType> propagate(const BoundParameters& start,
                               const propagator_options_t& options,
                               bool createFinalParameters = true) const;

  /// @brief Propagate track parameters to a target surface
  ///
  /// This function propagates the track parameters until the actor of the
  /// options aborts or the maximum number of steps is reached. The navigator
  /// steers to the target surface, and the @c DefaultActor stops there. An
  /// actor of your own must stop there itself.
  ///
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result containing the propagation status and the
  ///         final track parameters
  template <typename propagator_options_t>
  Result<ResultType> propagate(const BoundParameters& start,
                               const Surface& target,
                               const propagator_options_t& options) const;

  /// @brief Propagate track parameters with a propagator state of the caller
  ///
  /// This function initializes the state, runs the propagation and builds the
  /// result. The actor state stays in @p state, so the caller can read it
  /// afterwards in `state.actorState`.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  ///
  /// @param [in,out] state The propagator state object, from @c makeState
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Surface to propagate to, or nullptr if the propagation
  ///             has no target surface
  /// @param [in] createFinalParameters Whether to produce parameters at the end of the propagation
  ///
  /// @return Propagation result containing the propagation status and the
  ///         final track parameters
  template <typename propagator_state_t>
  Result<ResultType> propagate(propagator_state_t& state,
                               const BoundParameters& start,
                               const Surface* target = nullptr,
                               bool createFinalParameters = true) const;

  /// @brief Builds the propagator state object
  ///
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] options Propagation options
  ///
  /// @return Propagator state object
  template <typename propagator_options_t>
  State<propagator_options_t> makeState(
      const propagator_options_t& options) const;

  /// @brief Initialize the propagator state
  ///
  /// This function initializes the propagator state for a new propagation. The
  /// target surface goes to the navigator, which steers to it.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  ///
  /// @param [in,out] state The propagator state object
  /// @param [in] start Initial track parameters to propagate
  /// @param [in] target Surface to propagate to, or nullptr if the propagation
  ///             has no target surface
  ///
  /// @return Indication if the initialization was successful
  template <typename propagator_state_t>
  [[nodiscard]] Result<void> initialize(propagator_state_t& state,
                                        const BoundParameters& start,
                                        const Surface* target) const;

  /// @brief Run the propagation on an initialized state
  ///
  /// This function runs the propagation until the actor aborts or the
  /// maximum number of steps is reached.
  ///
  /// @tparam propagator_state_t Type of the propagator state with options
  ///
  /// @param [in,out] state the propagator state object
  ///
  /// @return Propagation result
  template <typename propagator_state_t>
  Result<void> propagate(propagator_state_t& state) const;

  /// @brief Builds the propagator result object
  ///
  /// This function creates the propagator result object from the propagator
  /// state object. The `result` is passed to pipe a potential error from the
  /// propagation call. The `createFinalParameters` flag is used to determine
  /// if the result should contain final track parameters.
  ///
  /// The actor state stays in @p state.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  ///
  /// @param [in,out] state Propagator state object
  /// @param [in] result Result of the propagation
  /// @param [in] createFinalParameters Whether to produce parameters at the end of the propagation
  /// @param [in] target Surface to be used for the final parameters, if createFinalParameters is true
  ///             If nullptr, the current surface of the navigator will be used
  ///
  /// @return Propagation result
  template <typename propagator_state_t>
  Result<ResultType> makeResult(propagator_state_t& state, Result<void> result,
                                bool createFinalParameters,
                                const Surface* target = nullptr) const;

 private:
  /// Refresh position, navigation direction, and geometric curvature from the
  /// current stepper state. Propagates magnetic-field lookup failures.
  template <typename propagator_state_t>
  Result<void> updateNavigationState(propagator_state_t& state) const;

  /// Implementation of propagation algorithm
  Stepper m_stepper;

  /// Implementation of navigator
  Navigator m_navigator;

  std::shared_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "Acts/Propagator/Propagator.ipp"
