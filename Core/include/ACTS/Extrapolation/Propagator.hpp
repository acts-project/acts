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
#include "ACTS/Extrapolation/ObserverList.hpp"
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
  struct Result
  {
    Result(const TrackParameters& startParameters,
           const Status&          s = Status::pUNSET)
      : endParameters(startParameters), status(s), m_tResults()
    {
    }

    TrackParameters endParameters;
    Status          status = Status::pUNSET;

    template <typename R>
    const R&
    get() const
    {
      return std::get<R>(m_tResults);
    }

    template <typename R>
    R&
    get()
    {
      return std::get<R>(m_tResults);
    }

  private:
    std::tuple<ExResult...> m_tResults;
  };

private:
  template <typename TrackParameters, typename... observers>
  struct observer_list_type_helper
  {
    template <typename... args>
    using result_type = Result<TrackParameters, args...>;

    typedef ObserverList<result_type, observers...> type;
  };

public:
  template <typename TrackParameters, typename... observers>
  using observer_list_t =
      typename observer_list_type_helper<TrackParameters, observers...>::type;

  /// @brief propagate track parameters
  template <typename TrackParameters, typename ObserverList>
  typename ObserverList::result_type
  propagate(const TrackParameters& start, const ObserverList& obsList)
  {
    typename ObserverList::result_type r(start, Status::pINPROGRESS);

    TrackParameters previous = start;
    TrackParameters current  = start;
    for (unsigned int i = 0; i < 1000; ++i) {
      current = m_impl.doStep(previous, 1 * units::_cm);
      obsList(current, previous, r);
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
