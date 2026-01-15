// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/Propagation.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/PropagatorResult.hpp"
#include "Acts/Propagator/PropagatorState.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Propagator/detail/ParameterTraits.hpp"
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
template <typename stepper_t, typename navigator_t = VoidNavigator>
class Propagator final
    : public std::conditional_t<
          SupportsBoundParameters_v<stepper_t>,
          detail::BasePropagatorHelper<Propagator<stepper_t, navigator_t>>,
          detail::PropagatorStub> {
  /// Re-define bound track parameters dependent on the stepper
  using StepperBoundTrackParameters =
      detail::stepper_bound_parameters_type_t<stepper_t>;
  static_assert(BoundTrackParametersConcept<StepperBoundTrackParameters>,
                "Stepper bound track parameters do not fulfill bound "
                "parameters concept.");

  using Jacobian = BoundMatrix;
  using BoundState = std::tuple<StepperBoundTrackParameters, Jacobian, double>;

  static_assert(StepperStateConcept<typename stepper_t::State>,
                "Stepper does not fulfill stepper concept.");
  static_assert(StepperConcept<stepper_t>,
                "Stepper does not fulfill stepper concept.");

  /// @brief Helper struct determining the state's type
  ///
  /// @tparam propagator_options_t Propagator options type
  ///
  /// This helper struct provides type definitions to extract the correct
  /// propagation state type from a given TrackParameter type and an
  /// ActorList.
  ///
  template <typename propagator_options_t>
  struct state_type_helper {
    using actor_list_t = typename propagator_options_t::actor_list_type;

    /// @brief Propagation state type for an arbitrary list of additional
    ///        propagation states
    ///
    /// @tparam args Parameter pack specifying additional propagation states
    ///
    template <typename... args>
    using this_state_type =
        PropagatorState<propagator_options_t, typename stepper_t::State,
                        typename navigator_t::State, args...>;

    /// @brief Propagation result type derived from a given action list
    using type = typename actor_list_t::template result_type<this_state_type>;
  };

  /// @brief Helper struct determining the result's type
  ///
  /// @tparam propagator_options_t Propagator options type
  ///
  /// This helper struct provides type definitions to extract the correct
  /// propagation result type from a given TrackParameter type and an
  /// ActorList.
  ///
  template <typename propagator_options_t>
  struct result_type_helper {
    using parameters_t = StepperBoundTrackParameters;
    using actor_list_t = typename propagator_options_t::actor_list_type;

    /// @brief Propagation result type for an arbitrary list of additional
    ///        propagation results
    ///
    /// @tparam args Parameter pack specifying additional propagation results
    ///
    template <typename... args>
    using this_result_type = PropagatorResult<parameters_t, args...>;

    /// @brief Propagation result type derived from a given action list
    using type = typename actor_list_t::template result_type<this_result_type>;
  };

 public:
  /// Type of the stepper in use
  using Stepper = stepper_t;

  /// Type of the navigator in use
  using Navigator = navigator_t;

  /// Type of the stepper options
  using StepperOptions = typename Stepper::Options;

  /// Type of the navigator options
  using NavigatorOptions = typename Navigator::Options;

  /// Type of the propagator options with actor list
  template <typename actor_list_t = ActorList<>>
  using Options =
      PropagatorOptions<StepperOptions, NavigatorOptions, actor_list_t>;

  /// Type of the stepper state
  using StepperState = typename Stepper::State;

  /// Type of the navigator state
  using NavigatorState = typename Navigator::State;

  /// Type of the propagation state derived from the propagator options
  /// @tparam propagator_options_t Type of the propagator options
  template <typename propagator_options_t>
  using State = typename state_type_helper<propagator_options_t>::type;

  /// Type of the propagation result derived from the propagator options
  /// @tparam propagator_options_t Type of the propagator options
  template <typename propagator_options_t>
  using ResultType = typename result_type_helper<propagator_options_t>::type;

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
  /// @param [in] createFinalParameters Whether to produce parameters at the end of the propagation
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t, typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  Result<ResultType<propagator_options_t>> propagate(
      const parameters_t& start, const propagator_options_t& options,
      bool createFinalParameters = true) const;

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
  Result<ResultType<propagator_options_t>> propagate(
      const parameters_t& start, const Surface& target,
      const propagator_options_t& options) const;

  /// @brief Builds the propagator state object
  ///
  /// This function creates the propagator state object from the initial track
  /// parameters and the propagation options.
  ///
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] options Propagation options
  ///
  /// @return Propagator state object
  template <typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  auto makeState(const propagator_options_t& options) const;

  /// @brief Builds the propagator state object
  ///
  /// This function creates the propagator state object from the initial track
  /// parameters, the target surface, and the propagation options.
  ///
  /// @tparam propagator_options_t Type of the propagator options
  /// @tparam target_aborter_t The target aborter type to be added
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagator state object
  template <typename propagator_options_t,
            typename target_aborter_t = SurfaceReached,
            typename path_aborter_t = PathLimitReached>
  auto makeState(const Surface& target,
                 const propagator_options_t& options) const;

  /// @brief Initialize the propagator state
  ///
  /// This function initializes the propagator state for a new propagation.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  /// @tparam path_aborter_t The path aborter type to be added
  ///
  /// @param [in,out] state The propagator state object
  /// @param [in] start Initial track parameters to propagate
  ///
  /// @return Indication if the initialization was successful
  template <typename propagator_state_t, typename parameters_t,
            typename path_aborter_t = PathLimitReached>
  [[nodiscard]] Result<void> initialize(propagator_state_t& state,
                                        const parameters_t& start) const;

  /// @brief Propagate track parameters
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @note Does not (yet) convert into the return_type of the propagation
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
  /// propagation call. The `options` are used to determine the type of the
  /// result object. The `createFinalParameters` flag is used to determine if
  /// the result should contain final track parameters.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] state Propagator state object
  /// @param [in] result Result of the propagation
  /// @param [in] options Propagation options
  /// @param [in] createFinalParameters Whether to produce parameters at the end of the propagation
  ///
  /// @return Propagation result
  template <typename propagator_state_t, typename propagator_options_t>
  Result<ResultType<propagator_options_t>> makeResult(
      propagator_state_t state, Result<void> result,
      const propagator_options_t& options, bool createFinalParameters) const;

  /// @brief Builds the propagator result object
  ///
  /// This function creates the propagator result object from the propagator
  /// state object. The `result` is passed to pipe a potential error from the
  /// propagation call. The `options` are used to determine the type of the
  /// result object.
  ///
  /// @tparam propagator_state_t Type of the propagator state object
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] state Propagator state object
  /// @param [in] result Result of the propagation
  /// @param [in] target Target surface of to propagate to
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result
  template <typename propagator_state_t, typename propagator_options_t>
  Result<ResultType<propagator_options_t>> makeResult(
      propagator_state_t state, Result<void> result, const Surface& target,
      const propagator_options_t& options) const;

  template <typename propagator_state_t>
  Result<NavigationTarget> getNextTarget(propagator_state_t& state) const;

  template <typename propagator_state_t>
  Result<void> performStep(propagator_state_t& state) const;

  template <typename propagator_state_t>
  Result<void> reachNextSurface(propagator_state_t& state) const;

  template <typename propagator_state_t>
  Result<bool> prePropagation(propagator_state_t& state) const;

  template <typename propagator_state_t>
  Result<void> postPropagation(propagator_state_t& state) const;

  template <typename propagator_state_t>
  class StepByStepPropagation {
   public:
    using StepByStepPropagationType = StepByStepPropagation<propagator_state_t>;
    using PropagatorType = Propagator<Stepper, Navigator>;
    using PropagationType = Propagation<PropagatorType, propagator_state_t>;

    class Iterator {
     public:
      Iterator() = default;
      Iterator(StepByStepPropagationType& propagation, bool done)
          : m_propagation(&propagation), m_done(done) {}
      Iterator(const Iterator&) = delete;
      Iterator(Iterator&&) = default;
      Iterator& operator=(const Iterator&) = delete;
      Iterator& operator=(Iterator&&) = default;

      bool operator==(const Iterator& other) const {
        return m_propagation == other.m_propagation && m_done == other.m_done;
      }

      std::tuple<const Result<void>&, PropagationType> operator*() const {
        return {m_lastResult, propagation()};
      }

      Iterator& operator++() {
        if (m_done || !m_lastResult.ok()) {
          return *this;
        }

        const Result<void> resStep = propagation().performStep();
        m_lastResult = resStep;
        m_done = propagation().state().terminatedNormally;

        return *this;
      }

     private:
      StepByStepPropagationType* m_propagation{nullptr};
      Result<void> m_lastResult = Result<void>::success();
      bool m_done{false};

      PropagationType propagation() const {
        return PropagationType(*m_propagation->m_propagator,
                               m_propagation->m_state);
      }
    };

    StepByStepPropagation(const PropagatorType& propagator,
                          propagator_state_t state)
        : m_propagator(&propagator), m_state(std::move(state)) {}

    Iterator begin() { return Iterator(*this, false); }
    Iterator end() { return Iterator(*this, true); }

   private:
    const PropagatorType* m_propagator{nullptr};
    propagator_state_t m_state;
  };

  template <typename parameters_t, typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  auto propagateStepByStep(const parameters_t& start,
                           const propagator_options_t& options) const {
    auto state = makeState<propagator_options_t, path_aborter_t>(options);
    using StateType = decltype(state);

    auto initRes =
        initialize<StateType, parameters_t, path_aborter_t>(state, start);
    // TODO check result

    auto preRes = prePropagation(state);
    // TODO check result

    return StepByStepPropagation<StateType>(*this, std::move(state));
  }

  template <typename propagator_state_t>
  class SurfaceBySurfacePropagation {
   public:
    using SurfaceBySurfacePropagationType =
        SurfaceBySurfacePropagation<propagator_state_t>;
    using PropagatorType = Propagator<Stepper, Navigator>;
    using PropagationType = Propagation<PropagatorType, propagator_state_t>;

    class Iterator {
     public:
      Iterator() = default;
      Iterator(SurfaceBySurfacePropagationType& propagation, bool done)
          : m_propagation(&propagation), m_done(done) {}
      Iterator(const Iterator&) = delete;
      Iterator(Iterator&&) = default;
      Iterator& operator=(const Iterator&) = delete;
      Iterator& operator=(Iterator&&) = default;

      bool operator==(const Iterator& other) const {
        return m_propagation == other.m_propagation && m_done == other.m_done;
      }

      std::tuple<const Result<void>&, const Surface*, PropagationType>
      operator*() const {
        return {m_lastResult, propagation().currentSurface(), propagation()};
      }

      Iterator& operator++() {
        if (m_done || !m_lastResult.ok()) {
          return *this;
        }

        const Result<void> resNextSurface = propagation().reachNextSurface();
        m_lastResult = resNextSurface.error();
        m_done = propagation().state().terminatedNormally;

        return *this;
      }

     private:
      SurfaceBySurfacePropagationType* m_propagation{nullptr};
      Result<void> m_lastResult = Result<void>::success();
      bool m_done{false};

      PropagationType propagation() const {
        return PropagationType(*m_propagation->m_propagator,
                               m_propagation->m_state);
      }
    };

    SurfaceBySurfacePropagation(const PropagatorType& propagator,
                                propagator_state_t state)
        : m_propagator(&propagator), m_state(std::move(state)) {}

    Iterator begin() { return Iterator(*this, false); }
    Iterator end() { return Iterator(*this, true); }

   private:
    const PropagatorType* m_propagator{nullptr};
    propagator_state_t m_state;
  };

  template <typename parameters_t, typename propagator_options_t,
            typename path_aborter_t = PathLimitReached>
  auto propagateSurfaceBySurface(const parameters_t& start,
                                 const propagator_options_t& options) const {
    auto state = makeState<propagator_options_t, path_aborter_t>(options);
    using StateType = decltype(state);

    auto initRes =
        initialize<StateType, parameters_t, path_aborter_t>(state, start);
    // TODO check result

    auto preRes = prePropagation(state);
    // TODO check result

    return SurfaceBySurfacePropagation<StateType>(*this, std::move(state));
  }

 private:
  /// Implementation of propagation algorithm
  Stepper m_stepper;

  /// Implementation of navigator
  Navigator m_navigator;

  std::shared_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }

  template <typename propagator_state_t, typename propagator_result_t>
  void moveStateToResult(propagator_state_t& state,
                         propagator_result_t& result) const;
};

}  // namespace Acts

#include "Acts/Propagator/Propagator.ipp"
