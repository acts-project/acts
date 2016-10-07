// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATION_PROPAGATOR_H
#define ACTS_EXTRAPOLATION_PROPAGATOR_H 1

#include <memory>
#include <type_traits>
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/detail/Extendable.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

namespace propagation {

  /// propagation direction relative to momentum
  enum Direction : int { backward = -1, forward = 1 };

  /// result status of track parameter propagation
  enum struct Status { SUCCESS, FAILURE, UNSET, INPROGRESS, WRONG_DIRECTION };

  /// result of propagation call
  template <typename TrackParameters, typename... ExResult>
  struct Result : private detail::Extendable<ExResult...>
  {
    using detail::Extendable<ExResult...>::get;

    std::unique_ptr<TrackParameters> endParameters = nullptr;
    Status                           status        = Status::UNSET;
    unsigned int                     steps         = 0;
    double                           pathLength    = 0.;

    operator bool() const
    {
      return (endParameters && status == Status::SUCCESS);
    }
  };

  /// @brief propagator for particles in a magnetic field
  ///
  /// @tparam Impl implementation of the propagation algorithm
  template <typename Impl>
  class Propagator final
  {
  public:
    /// @brief options for propagate() call
    template <typename ObserverList = ObserverList<>,
              typename AbortList    = AbortList<>>
    struct Options
    {
      /// propagation direction
      Direction direction = forward;

      /// maximum number of steps for one propagate() call
      unsigned int max_steps = 1000;

      /// minimum step size
      double min_step_size = 0.1 * units::_mm;

      /// maximum step size
      double max_step_size = 1 * units::_m;

      /// maximum step size
      double max_path_length = 5 * units::_m;

      /// list of observers
      ObserverList observer_list;

      /// list of stop conditions
      AbortList stop_conditions;
    };

    /// @brief copy implementation object
    template <typename T = Impl>
    explicit Propagator(
        typename std::enable_if_t<std::is_copy_constructible<T>::value,
                                  const T&> impl)
      : m_impl(impl)
    {
    }

    /// @brief move implementation object
    template <typename T = Impl>
    explicit Propagator(
        typename std::enable_if_t<std::is_move_constructible<T>::value, T&&>
            impl
        = Impl())
      : m_impl(std::move(impl))
    {
    }

    /// @brief copy constructor
    template <typename T = Impl>
    Propagator(typename std::enable_if_t<std::is_copy_constructible<T>::value,
                                         const Propagator<T>&> copy)
      : m_impl(copy.m_impl)
    {
    }

    /// @brief move constructor
    template <typename T = Impl>
    Propagator(typename std::enable_if_t<std::is_move_constructible<T>::value,
                                         Propagator<T>&&> rhs)
      : m_impl(std::move(rhs.m_impl))
    {
    }

    /// @brief copy assignment operator
    template <typename T = Impl>
    typename std::enable_if_t<std::is_copy_constructible<T>::value,
                              const Propagator<T>&>
    operator=(const Propagator<T>& copy)
    {
      m_impl = copy.m_impl;
      return *this;
    }

    /// @brief move assignment operator
    template <typename T = Impl>
    typename std::enable_if_t<std::is_move_constructible<T>::value,
                              const Propagator<T>&>
    operator=(Propagator<T>&& copy)
    {
      m_impl = std::move(copy.m_impl);
      return *this;
    }

  private:
    template <typename TrackParameters, typename ObserverList>
    struct result_type_helper
    {
      template <typename... args>
      using this_result_type = Result<TrackParameters, args...>;

      typedef typename ObserverList::template result_type<this_result_type>
          type;
    };

    template <typename T, typename ObsList>
    using obs_list_result_t = typename result_type_helper<T, ObsList>::type;

  public:
    /// @brief propagate track parameters
    template <typename TrackParameters,
              typename ObserverList,
              typename AbortList>
    obs_list_result_t<
        typename Impl::template return_parameter_type<TrackParameters>,
        ObserverList>
    propagate(const TrackParameters& start,
              const Options<ObserverList, AbortList>& options) const
    {
      typedef typename Impl::template cache_type<TrackParameters> cache_type;
      typedef typename Impl::template return_parameter_type<TrackParameters>
          return_parameter_type;
      typedef obs_list_result_t<return_parameter_type, ObserverList>
          result_type;

      static_assert(std::is_copy_constructible<return_parameter_type>::value,
                    "return track parameter type must be copy-constructible");

      result_type  r;
      const double signed_pathLimit
          = options.direction * options.max_path_length;
      double                stepMax = options.direction * options.max_step_size;
      cache_type            propagation_cache(start);
      return_parameter_type previous = m_impl.convert(propagation_cache);
      for (; r.steps < options.max_steps; ++r.steps) {
        r.pathLength += m_impl.step(propagation_cache, stepMax);
        return_parameter_type current = m_impl.convert(propagation_cache);
        options.observer_list(current, previous, r);
        if (fabs(r.pathLength) >= options.max_path_length
            || options.stop_conditions(r, current, stepMax)) {
          r.endParameters = std::make_unique<return_parameter_type>(
              m_impl.convert(propagation_cache));
          r.status = Status::SUCCESS;
          break;
        }

        if (fabs(stepMax) > fabs(signed_pathLimit - r.pathLength))
          stepMax = signed_pathLimit - r.pathLength;

        previous = current;
      }

      return r;
    }

    /// @brief propagate track parameters
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
      typedef typename Impl::template cache_type<TrackParameters, Surface>
          cache_type;
      typedef typename Impl::template step_parameter_type<TrackParameters>
          step_parameter_type;
      typedef typename Impl::template return_parameter_type<TrackParameters,
                                                            Surface>
          return_parameter_type;
      typedef obs_list_result_t<return_parameter_type, ObserverList>
          result_type;

      static_assert(std::is_copy_constructible<return_parameter_type>::value,
                    "return track parameter type must be copy-constructible");

      result_type  r;
      const double signed_pathLimit
          = options.direction * options.max_path_length;
      cache_type          cache(start);
      step_parameter_type previous = m_impl.convert(cache);
      double              distance
          = m_impl.distance(target, cache.position(), cache.direction());

      if (distance * options.direction < 0) {
        r.status = Status::WRONG_DIRECTION;
        return r;
      }

      double stepMax = options.direction * options.max_step_size;
      if (fabs(stepMax) > fabs(distance)) stepMax = distance;

      for (; r.steps < options.max_steps; ++r.steps) {
        r.pathLength += m_impl.step(cache, stepMax);
        step_parameter_type current = m_impl.convert(cache);
        options.observer_list(current, previous, r);
        distance = m_impl.distance(target, cache.position(), cache.direction());
        if (fabs(distance) < 1 * units::_um
            || fabs(r.pathLength) >= options.max_path_length
            || options.stop_conditions(r, current, stepMax)) {
          r.endParameters = std::make_unique<return_parameter_type>(
              m_impl.convert(cache, target));
          r.status = Status::SUCCESS;
          break;
        }

        if (fabs(stepMax) > fabs(signed_pathLimit - r.pathLength))
          stepMax = signed_pathLimit - r.pathLength;

        if (fabs(stepMax) > fabs(distance)) stepMax = distance;

        previous = current;
      }

      return r;
    }

  private:
    /// implementation of propagation algorithm
    Impl m_impl;
  };

}  // namespace propagation

}  // namespace Acts

#endif  // ACTS_EXTRAPOLATION_PROPAGATOR_H
