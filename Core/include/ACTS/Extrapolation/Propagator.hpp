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
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/detail/Extendable.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

/// @brief propagator for particles in a magnetic field
///
/// @tparam Impl implementation of the propagation algorithm
template <typename Impl>
class Propagator final
{
public:
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
  template <typename TrackParameters, typename ObserverList, typename AbortList>
  obs_list_result_t<TrackParameters, ObserverList>
  propagate(const TrackParameters& start,
            const ObserverList&    obsList,
            const AbortList&       conditions)
  {
    typedef obs_list_result_t<TrackParameters, ObserverList> result_type;
    result_type     r(start, Status::pINPROGRESS);
    double          stepMax  = 1 * units::_m;
    TrackParameters previous = start;
    TrackParameters current  = start;
    for (unsigned int i = 0; i < 1000; ++i) {
      current = m_impl.doStep(previous, stepMax);
      obsList(current, previous, r);
      if (conditions(current, r, stepMax)) break;
      previous = current;
    }

    r.endParameters = current;
    r.status        = Status::pSUCCESS;

    return r;
  }

private:
  /// implementation of propagation algorithm
  Impl m_impl;
};

}  // end of namespace

#endif  // ACTS_EXTRAPOLATION_PROPAGATOR_H
