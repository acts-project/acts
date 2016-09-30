// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTRAPOLATION_PROPAGATOR_H
#define ACTS_EXTRAPOLATION_PROPAGATOR_H 1

#include <type_traits>
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/Direction.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/detail/Extendable.hpp"
#include "ACTS/Extrapolation/detail/step_caller.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

/// @brief propagator for particles in a magnetic field
///
/// @tparam Impl implementation of the propagation algorithm
template <typename Impl>
class Propagator final
{
public:
  /// @brief options for propagate() call
  template <Direction             = forward,
            typename ObserverList = ObserverList<>,
            typename AbortList    = AbortList<>>
  struct Options
  {
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
      typename std::enable_if_t<std::is_copy_constructible<T>::value, const T&>
          impl)
    : m_impl(impl)
  {
  }

  /// @brief move implementation object
  template <typename T = Impl>
  explicit Propagator(
      typename std::enable_if_t<std::is_move_constructible<T>::value, T&&> impl
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

  enum struct Status { pSUCCESS, pFAILURE, pUNSET, pINPROGRESS };

  template <typename TrackParameters, typename... ExResult>
  struct Result : private detail::Extendable<ExResult...>
  {
    Result(const TrackParameters& startParameters,
           const Status&          s = Status::pUNSET)
      : detail::Extendable<ExResult...>()
      , endParameters(startParameters)
      , status(s)
    {
    }

    using detail::Extendable<ExResult...>::get;

    TrackParameters endParameters;
    Status          status = Status::pUNSET;
  };

private:
  template <typename TrackParameters, typename ObserverList>
  struct result_type_helper
  {
    template <typename... args>
    using this_result_type = Result<TrackParameters, args...>;

    typedef typename ObserverList::template result_type<this_result_type> type;
  };

  template <typename T, typename ObsList>
  using obs_list_result_t = typename result_type_helper<T, ObsList>::type;

public:
  /// @brief propagate track parameters
  template <typename TrackParameters,
            Direction direction,
            typename ObserverList,
            typename AbortList>
  obs_list_result_t<
      typename Impl::template return_parameter_type<TrackParameters>,
      ObserverList>
  propagate(const TrackParameters& start,
            const Options<direction, ObserverList, AbortList>& options) const
  {
    typedef typename Impl::template cache_type<TrackParameters> cache_type;
    typedef typename Impl::template return_parameter_type<TrackParameters>
        return_parameter_type;
    typedef obs_list_result_t<return_parameter_type, ObserverList> result_type;
    typedef detail::step_caller<Impl, cache_type, direction> step_caller;

    static_assert(std::is_copy_constructible<return_parameter_type>::value,
                  "return track parameter type must be copy-constructible");

    result_type           r(start, Status::pINPROGRESS);
    double                stepMax = options.max_step_size;
    cache_type            propagation_cache(start, direction);
    return_parameter_type previous
        = m_impl.convert(propagation_cache, direction);
    double pathLength = 0;
    for (unsigned int i = 0; i < options.max_steps; ++i) {
      pathLength += step_caller::step(m_impl, propagation_cache, stepMax);
      return_parameter_type current
          = m_impl.convert(propagation_cache, direction);
      options.observer_list(current, previous, r);
      if (pathLength >= options.max_path_length
          || options.stop_conditions(r, current, stepMax)) {
        r.endParameters = m_impl.convert(propagation_cache, direction);
        r.status        = Status::pSUCCESS;
        break;
      }

      if (stepMax > options.max_path_length - pathLength)
        stepMax = options.max_path_length - pathLength;

      previous = current;
    }

    return r;
  }

private:
  /// implementation of propagation algorithm
  Impl m_impl;
};

}  // end of namespace

#endif  // ACTS_EXTRAPOLATION_PROPAGATOR_H
