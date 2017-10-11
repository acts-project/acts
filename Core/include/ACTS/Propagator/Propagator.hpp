// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATION_PROPAGATOR_H
#define ACTS_EXTRAPOLATION_PROPAGATOR_H 1

#include <cmath>
#include <memory>
#include <type_traits>
#include "ACTS/Propagator/AbortList.hpp"
#include "ACTS/Propagator/ObserverList.hpp"
#include "ACTS/Propagator/detail/Extendable.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace propagation {

  /// Propagation direction, relative to momentum
  enum Direction : int { backward = -1, forward = 1 };

  /// Result status of track parameter propagation
  enum struct Status { SUCCESS, FAILURE, UNSET, INPROGRESS, WRONG_DIRECTION };

  /// @brief Simple class holding result of propagation call
  ///
  /// @tparam TrackParameters Type of final track parameters
  ///
  /// @tparam ExResult        Parameter pack for additional propagation
  ///                         quantities
  ///
  template <typename TrackParameters, typename... ExResult>
  struct Result : private detail::Extendable<ExResult...>
  {
    /// Constructor from initial propagation status
    Result(Status s = Status::UNSET)
      : detail::Extendable<ExResult...>(), status(s)
    {
    }

    /// Accessor to additional propagation quantities
    using detail::Extendable<ExResult...>::get;

    /// Final track parameters
    std::unique_ptr<const TrackParameters> endParameters = nullptr;

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
    operator bool() const
    {
      return (endParameters && status == Status::SUCCESS);
    }
  };

  /// @brief Propagator for particles in a magnetic field
  ///
  /// @tparam Impl Implementation of the propagation algorithm
  ///
  /// This Propagator class serves as high-level steering code for propagating
  /// track parameters. The actual implementation of the propagation has to be
  /// implemented in the Impl object, which has to provide the following:
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
  template <typename Impl>
  class Propagator final
  {
  public:
    
    // Type of cache object used by the propagation implementation
    typedef typename Impl::template cache_type<TrackParameters> cache_type;
    
    /// @brief Options for propagate() call
    ///
    /// @tparam ObserverList List of observer types called after each
    ///                      propagation step with the current propagation
    ///                      cache
    ///
    /// @tparam AbortList    List of abort conditions tested after each
    ///                      propagation step using the current propagation
    ///                      cache
    ///
    template <typename ObserverList = ObserverList<>,
              typename AbortList    = AbortList<>>
    struct Options
    {
      /// Propagation direction
      Direction direction = forward;

      /// Maximum number of steps for one propagate() call
      unsigned int max_steps = 1000;

      /// Required distance to surface 
      double target_surface_distance = 1 * units::_um;

      /// Absolute minimum step size
      double min_step_size = 0.1 * units::_mm;

      /// Absolute maximum step size
      double max_step_size = 1 * units::_m;

      /// Absolute maximum path length
      double max_path_length = 5 * units::_m;

      /// List of observers
      ObserverList observer_list;

      /// List of abort conditions
      AbortList stop_conditions;
    };

    /// Constructor from implementation object
    explicit Propagator(Impl impl) : m_impl(std::move(impl)) {}

  private:
    /// @brief Helper struct determining the result's type
    ///
    /// @tparam TrackParameters Type of final track parameters
    /// @tparam ObserverList    List of propagation observer types
    ///
    /// This helper struct provides type definitions to extract the correct
    /// propagation result type from a given TrackParameter type and an
    /// ObserverList.
    ///
    template <typename TrackParameters, typename ObserverList>
    struct result_type_helper
    {
      /// @brief Propagation result type for an arbitrary list of additional
      ///        propagation results
      ///
      /// @tparam args Parameter pack specifying additional propagation results
      ///
      template <typename... args>
      using this_result_type = Result<TrackParameters, args...>;

      /// @brief Propagation result type derived from a given observer list
      typedef typename ObserverList::template result_type<this_result_type>
          type;
    };

    /// @brief Short-hand type definition for propagation result derived from
    ///        an observer list
    ///
    /// @tparam T       Type of the final track parameters
    /// @tparam ObsList List of propagation observer types
    ///
    template <typename T, typename ObsList>
    using obs_list_result_t = typename result_type_helper<T, ObsList>::type;
    
  public:
        
    /// @brief Propagate track parameters - User method
    ///
    /// This function performs the propagation of the track parameters using the
    /// internal implementation object, until at least one abort condition is
    /// fulfilled or the maximum number of steps/path length provided in the
    /// propagation options is reached.
    ///
    /// @tparam TrackParameters Type of initial track parameters to propagate
    /// @tparam ObserverList    Type list of observers
    /// @tparam AbortList       Type list of abort conditions
    ///
    /// @param [in] start   Initial track parameters to propagate
    /// @param [in] options Propagation options
    ///
    /// @return Propagation result containing the propagation status, final
    ///         track parameters, and output of observers (if they produce any)
    ///
    template <typename TrackParameters,
              typename ObserverList,
              typename AbortList>
    obs_list_result_t<
        typename Impl::template return_parameter_type<TrackParameters>,
        ObserverList>
    propagate(const TrackParameters& start,
              const Options<ObserverList, AbortList>& options) const
    {

      // Type of track parameters produced by the propagation
      typedef typename Impl::template return_parameter_type<TrackParameters>
          return_parameter_type;

      // Type of the full propagation result, including output from observers
      typedef obs_list_result_t<return_parameter_type, ObserverList>
          result_type;

      static_assert(std::is_copy_constructible<return_parameter_type>::value,
                    "return track parameter type must be copy-constructible");

      // Initialize the propagation result object
      result_type  r(Status::INPROGRESS);

      // Compute the signed path limit and maximum step size
      const double signed_pathLimit
          = options.direction * options.max_path_length;
      double stepMax = options.direction * options.max_step_size;

      // Initialize the internal propagation cache
      cache_type cache(start);

      // Propagation loop
      for (; r.steps < options.max_steps; ++r.steps) {
        // Perform a propagation step (which can alter the step size)
        r.pathLength += m_impl.step(cache, stepMax);
        // Call the observers with the current and previous track parameters,
        // and let them fill in some propagation results
        options.observer_list(cache, r);
        // Is it time to stop the propagation?
        if (std::abs(r.pathLength) >= options.max_path_length
            || options.stop_conditions(r, cache, stepMax)) {
          // Compute the final results and mark the propagation as successful
          r.endParameters = std::make_unique<const return_parameter_type>(
              m_impl.convert(cache));
          r.status = Status::SUCCESS;
          break;
        }
        // Adjust the step size so that we cannot go above the path limit
        if (std::abs(stepMax) > std::abs(signed_pathLimit - r.pathLength))
          stepMax = signed_pathLimit - r.pathLength;
      }
      return r;      
    }

    /// @brief Propagate track parameters - Expert method with propagation cache
    ///
    /// This function performs the propagation of the track parameters according
    /// to the internal implementation object until at least one abort condition
    /// is fulfilled, the destination surface is hit or the maximum number of
    /// steps/path length as given in the propagation options is reached.
    ///
    /// It does check/re-use the propgation cache as much as possible for 
    /// performance reasons
    ///
    /// @tparam TrackParameters Type of initial track parameters to propagate
    /// @tparam Surface         Type of target surface
    /// @tparam ObserverList    Type list of observers
    /// @tparam AbortList       Type list of abort conditions
    ///
    /// @param [in] cache Stepper cache built/updated from the start parameters
    /// @param [in] target Target surface of to propagate to 
    /// @param [in] options Propagation options
    ///
    /// @return Propagation result containing the propagation status, final
    ///         track parameters, and output of observers (if they produce any)
    ///
    template <typename TrackParameters,
              typename Surface,
              typename ObserverList,
              typename AbortList>
    obs_list_result_t<
        typename Impl::template return_parameter_type<TrackParameters, Surface>,
        ObserverList>
    propagate_with_cache(cache_type& cache,
                         const TrackParameters& start, 
                         const Surface& target,
                         const Options<ObserverList, AbortList>& options) const
    {
      
      // Type of track parameters produced at the end of the propagation
      typedef typename Impl::template return_parameter_type<TrackParameters,
                                                            Surface>
          return_parameter_type;

      // Type of the full propagation result, including output from observers
      typedef obs_list_result_t<return_parameter_type, ObserverList>
          result_type;

      static_assert(std::is_copy_constructible<return_parameter_type>::value,
                    "return track parameter type must be copy-constructible");

      // Initialize the propagation result object
      result_type  r(Status::INPROGRESS);

      // Compute the signed path limit and maximum step size
      const double signed_pathLimit
          = options.direction * options.max_path_length;
      double stepMax = options.direction * options.max_step_size;

      // Compute the distance to the target surface
      double distance
          = m_impl.distance(target, cache.position(), cache.direction());

      // Check that we are propagating towards the surface
      if (distance * options.direction < 0) {
        r.status = Status::WRONG_DIRECTION;
        return r;
      }

      // Avoid going beyond the target surface on the first step
      if (std::abs(stepMax) > std::abs(distance)) stepMax = distance;

      // Propagation loop
      for (; r.steps < options.max_steps; ++r.steps) {
        // Perform a propagation step (which can alter the step size)
        r.pathLength += m_impl.step(cache, stepMax);

        // Call the observers with the current and previous track parameters,
        // and let them fill in some propagation results
        options.observer_list(cache, r);

        // Update the distance to the target surface
        distance = m_impl.distance(target, cache.position(), cache.direction());

        // Is it time to stop the propagation ?
        if (std::abs(distance) < options.target_surface_distance 
            || std::abs(r.pathLength) >= options.max_path_length
            || options.stop_conditions(r, cache, stepMax)) {
          // Compute the final results and mark the propagation as successful
          r.endParameters = std::make_unique<const return_parameter_type>(
              m_impl.convert(cache, target));
          r.status = Status::SUCCESS;
          break;
        }
        // Adjust the step size so that we cannot go above the path limit
        if (std::abs(stepMax) > std::abs(signed_pathLimit - r.pathLength))
          stepMax = signed_pathLimit - r.pathLength;
        // Adjust the step size so that we cannot cross the target surface
        if (std::abs(stepMax) > std::abs(distance)) stepMax = distance;
      }
      return r;
    }

    /// @brief Propagate track parameters - User method
    ///
    /// This function performs the propagation of the track parameters according
    /// to the internal implementation object until at least one abort condition
    /// is fulfilled, the destination surface is hit or the maximum number of
    /// steps/path length as given in the propagation options is reached.
    ///
    /// A stepper cache object is built internally for this call and the
    /// Expert method with the cache call signature is called.
    ///
    /// @tparam TrackParameters Type of initial track parameters to propagate
    /// @tparam Surface         Type of target surface
    /// @tparam ObserverList    Type list of observers
    /// @tparam AbortList       Type list of abort conditions
    ///
    /// @param [in] start Initial track parameters to propagate
    /// @param [in] target Target surface of to propagate to 
    /// @param [in] options Propagation options
    ///
    /// @return Propagation result containing the propagation status, final
    ///         track parameters, and output of observers (if they produce any)
    ///
    template <typename TrackParameters,
              typename Surface,
              typename ObserverList,
              typename AbortList>
    obs_list_result_t<
        typename Impl::template return_parameter_type<TrackParameters, Surface>,
        ObserverList>
    propagate(const TrackParameters& start,
              const Surface&         target,
              const Options<ObserverList, AbortList>& options) const
    {
      // Initialize the internal propagation cache
      cache_type  cache(start);
      return propagate_with_cache(cache, start, target, options);
    }

  private:
    /// implementation of propagation algorithm
    Impl m_impl;
  };

}  // namespace propagation

}  // namespace Acts

#endif  // ACTS_EXTRAPOLATION_PROPAGATOR_H
