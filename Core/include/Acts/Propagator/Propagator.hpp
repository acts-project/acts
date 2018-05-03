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
#include "Acts/Propagator/detail/standard_abort_conditions.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// Result status of track parameter propagation
enum struct Status { SUCCESS, FAILURE, UNSET, IN_PROGRESS, WRONG_DIRECTION };

/// @brief Simple class holding result of propagation call
///
/// @tparam parameters_t Type of final track parameters
///
/// @tparam result_list  Result pack for additional propagation
///                      quantities
///
template <typename parameters_t, typename... result_list>
struct Result : private detail::Extendable<result_list...>
{
  /// Constructor from initial propagation status
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

/// @brief Propagator for particles (optionally in a magnetic field)
///
/// The Propagator works with two cache objects:
///  - a propgator cache for object navigation and screen output
///  - a stepper cache for the actual transport caching (pos,dir,field)
///
/// @tparam stepper_t stepper imentation of the propagation algorithm
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
///   cache object
/// - a type mapping for: (initial track parameter type and destination
///   surface type) -> type of internal cache object
///
template <typename stepper_t>
class Propagator final
{
public:
  /// Type of cache object used by the propagation implementation
  typedef typename stepper_t::template cache_type<TrackParameters> StepperCache;

  /// Typedef the PathLimitReached aborter
  typedef detail::PathLimitReached PathLimitReached;

  /// @brief Options for propagate() call
  ///
  /// @tparam action_list_t List of action types called after each
  ///                   propagation step with the current propagation
  ///                   cache
  ///
  /// @tparam aborter_list_t  List of abort conditions tested after each
  ///                   propagation step using the current propagation
  ///                   cache
  ///
  template <typename action_list_t  = ActionList<>,
            typename aborter_list_t = AbortList<>>
  struct Options
  {

    /// Propagation direction
    NavigationDirection direction = forward;

    /// Maximum number of steps for one propagate() call
    unsigned int maxSteps = 1000;

    /// Required tolerance to reach target (surface, pathlength)
    double targetTolerance = s_onSurfaceTolerance;

    /// Absolute maximum step size
    double maxStepSize = 1 * units::_m;

    /// Absolute maximum path length
    double maxPathLength = std::numeric_limits<double>::max();

    /// Debugging option, this steers the output stream
    bool debug = false;

    /// List of actions
    action_list_t actionList;

    /// List of abort conditions
    aborter_list_t stopConditions;
  };

  /// Constructor from implementation object
  explicit Propagator(stepper_t stepper) : m_stepper(std::move(stepper)) {}

private:
  /// @brief private Propagator cache for navigation and debugging
  ///
  /// This struct holds the common cache information for propagating
  /// which is independent of the actual stepper implementation.
  struct PropagatorCache
  {

    /// Navigation cache: the start surface
    const Surface* startSurface = nullptr;

    /// Navigation cache: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation cache: the target surface
    const Surface* targetSurface = nullptr;
    bool           targetReached = false;

    /// Navigation cache : a break has been detected
    bool navigationBreak = false;

    /// Debug output steering
    /// - the string where debug messages are stored (optionally)
    /// - it also has some formatting options
    bool        debug         = false;
    std::string debugString   = "";
    size_t      debugPfxWidth = 30;
    size_t      debugMsgWidth = 50;
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
    typedef typename action_list_t::template result_type<this_result_type> type;
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

  /// @brief Propagate track parameters - Private method with cache
  ///
  /// This function performs the propagation of the track parameters according
  /// to the internal implementation object until at least one abort condition
  /// is fulfilled, the destination surface is hit or the maximum number of
  /// steps/path length as given in the propagation options is reached.
  ///
  /// @note Does not (yet) convert into  the return_type of the propagation
  ///
  /// @tparam result_t Type of the result object for this propagation
  /// @tparam action_list_t  Type list of actions, type ActionList<>
  /// @tparam aborter_list_t Type list of abort conditions, type AbortList<>
  /// @tparam internal_aborter_list_t additional internal aborters
  ///
  /// @param [in,out] Result of the propagation
  /// @param [in,out] Cache Stepper cache built/updated from the start
  /// parameters
  /// @param [in] Target Target surface of to propagate to
  /// @param [in] Options Propagation options
  ///
  /// @return Propagation Status
  template <typename result_t,
            typename action_list_t,
            typename aborter_list_t,
            typename internal_aborter_list_t>
  Status
  propagate_(result_t&        result,
             PropagatorCache& pCache,
             StepperCache&    sCache,
             const Options<action_list_t, aborter_list_t>& options,
             const internal_aborter_list_t& interalAborters) const
  {

    // Pre-stepping call to the abort list
    debugLog(options.debug, pCache, [&] {
      return std::string("Calling pre-stepping aborters.");
    });
    if (interalAborters(result, pCache, sCache)) return Status::FAILURE;

    // Pre-stepping call to the action list
    debugLog(options.debug, pCache, [&] {
      return std::string("Calling pre-stepping action list.");
    });
    options.actionList(pCache, sCache, result);

    // Propagation loop : stepping
    for (; result.steps < options.maxSteps; ++result.steps) {
      // Perform a propagation step - it only takes the stepping cache
      result.pathLength += m_stepper.step(sCache);
      // Call the actions, can (& will likely) modify cache
      debugLog(options.debug, pCache, [&] {
        return std::string("Calling action list on single step.");
      });
      options.actionList(pCache, sCache, result);
      // Call the stop_conditions and the internal stop conditions
      // break condition triggered, but still count the step
      debugLog(options.debug, pCache, [&] {
        return std::string("Calling aborters on single step.");
      });
      if (options.stopConditions(result, pCache, sCache)
          || interalAborters(result, pCache, sCache)) {
        ++result.steps;
        break;
      }
    }

    // Post-stepping call to the action list
    debugLog(options.debug, pCache, [&] {
      return std::string("Calling post-stepping action list.");
    });
    options.actionList(pCache, sCache, result);
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
  /// @tparam action_list_t       Type list of actions, type ActionList<>
  /// @tparam aborter_list_t        Type list of abort conditions, type
  /// AbortList<>
  ///
  /// @param [in] start   Initial track parameters to propagate
  /// @param [in] options Propagation options
  ///
  /// @return Propagation result containing the propagation status, final
  ///         track parameters, and output of actions (if they produce any)
  ///
  template <typename parameters_t,
            typename action_list_t,
            typename aborter_list_t>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t>,
      action_list_t>
  propagate(const parameters_t& start,
            const Options<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced by the propagation
    typedef typename stepper_t::template return_parameter_type<parameters_t>
        return_parameter_type;

    // Type of the full propagation result, including output from actions
    typedef action_list_t_result_t<return_parameter_type, action_list_t> Result;

    static_assert(std::is_copy_constructible<return_parameter_type>::value,
                  "return track parameter type must be copy-constructible");

    // Get the reference surface for navigation
    const auto& startSurface = start.referenceSurface();

    // Initialize the propagation result object
    Result result(Status::IN_PROGRESS);

    // Initialize the interal propagator cache
    PropagatorCache pCache;
    pCache.startSurface = &startSurface;
    pCache.debug        = options.debug;

    // Initialize the internal stepper cache
    StepperCache sCache(start, options.direction, options.maxStepSize);

    // Internal Abort list
    AbortList<PathLimitReached> interalAborters;
    // configure the aborter
    auto& pathLimitAbort = interalAborters.template get<PathLimitReached>();
    pathLimitAbort.signedPathLimit
        = std::abs(options.maxPathLength) * options.direction;
    pathLimitAbort.tolerance = options.targetTolerance;
    pathLimitAbort.debug     = options.debug;

    // Perform the actual propagation & check it's outcome
    if (propagate_(result, pCache, sCache, options, interalAborters)
        != Status::IN_PROGRESS) {
      debugLog(options.debug, pCache, [&] {
        return std::string("Propagation was not successful.");
      });
    } else {
      /// Convert into the return type
      result.endParameters = std::make_unique<const return_parameter_type>(
          m_stepper.convert(sCache));
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
  /// @tparam Surface         Type of target surface
  /// @tparam action_list_t       Type list of actions
  /// @tparam aborter_list_t        Type list of abort conditions
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
            typename aborter_list_t>
  action_list_t_result_t<
      typename stepper_t::template return_parameter_type<parameters_t,
                                                         surface_t>,
      action_list_t>
  propagate(const parameters_t& start,
            const surface_t&    target,
            const Options<action_list_t, aborter_list_t>& options) const
  {

    // Type of track parameters produced at the end of the propagation
    typedef typename stepper_t::template return_parameter_type<parameters_t,
                                                               surface_t>
        return_parameter_type;

    // Get the reference surface for navigation
    const auto& startSurface = start.referenceSurface();

    // Initialize the interal propagator cache
    PropagatorCache pCache;
    pCache.startSurface = &startSurface;
    pCache.debug        = options.debug;

    // Initialize the internal stepper cache
    StepperCache sCache(start, options.direction, options.maxStepSize);

    // Type of the full propagation result, including output from actions
    typedef action_list_t_result_t<return_parameter_type, action_list_t> Result;

    // Initialize the propagation result object
    Result result(Status::IN_PROGRESS);

    static_assert(std::is_copy_constructible<return_parameter_type>::value,
                  "return track parameter type must be copy-constructible");

    // Target surface abort condition with tolerance
    typedef detail::SurfaceReached<surface_t> targetReached;

    // Internal Abort list
    AbortList<targetReached, PathLimitReached> interalAborters;
    // configure the aborters
    auto& target_abort     = interalAborters.template get<targetReached>();
    target_abort.surface   = &target;
    target_abort.direction = options.direction;
    target_abort.tolerance = options.targetTolerance;
    target_abort.debug     = options.debug;

    auto& pathLimitAbort = interalAborters.template get<PathLimitReached>();
    pathLimitAbort.signedPathLimit
        = std::abs(options.maxPathLength) * options.direction;
    pathLimitAbort.tolerance = options.targetTolerance;
    pathLimitAbort.debug     = options.debug;

    // Perform the actual propagation
    if (propagate_(result, pCache, sCache, options, interalAborters)
        != Status::IN_PROGRESS) {
      debugLog(options.debug, pCache, [&] {
        return std::string("Propagation was not successful.");
      });
    } else {
      // Compute the final results and mark the propagation as successful
      result.endParameters = std::make_unique<const return_parameter_type>(
          m_stepper.convert(sCache, target));
      result.status = Status::SUCCESS;
    }
    return result;
  }

private:
  /// implementation of propagation algorithm
  stepper_t m_stepper;

  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the cache.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @param cache the stepper cache for the debug flag, prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_cache_t>
  void
  debugLog(bool                         debug,
           propagator_cache_t&          pCache,
           std::function<std::string()> logAction) const
  {
    if (debug) {
      std::stringstream dstream;
      dstream << "|->" << std::setw(pCache.debugPfxWidth);
      dstream << "Propagator"
              << " | ";
      dstream << std::setw(pCache.debugMsgWidth) << logAction() << '\n';
      pCache.debugString += dstream.str();
      std::cout << dstream.str();
    }
  }
};

}  // namespace Acts
<<<<<<< HEAD:Core/include/Acts/Propagator/Propagator.hpp
=======

#endif  // ACTS_EXTRAPOLATION_PROPAGATOR_H
>>>>>>> 167ab66c... clang-format after message-format:Core/include/ACTS/Propagator/Propagator.hpp
