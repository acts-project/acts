// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ABORT_LIST_HPP
#define ACTS_ABORT_LIST_HPP 1

#include "ACTS/Propagator/ObserverList.hpp"
#include "ACTS/Propagator/detail/Extendable.hpp"
#include "ACTS/Propagator/detail/abort_condition_signature_check.hpp"
#include "ACTS/Propagator/detail/abort_list_implementation.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"

namespace Acts {

template <typename... conditions>
struct AbortList : private detail::Extendable<conditions...>
{
private:
  static_assert(not detail::has_duplicates_v<conditions...>,
                "same abort conditions type specified several times");

  // clang-format off
  typedef detail::type_collector_t<detail::observer_type_extractor, conditions...> observers;
  // clang-format on

  using detail::Extendable<conditions...>::tuple;

public:
  /// default constructor
  AbortList() = default;

  /// constructor of the Extendable
  AbortList(const std::tuple<conditions...>& extensions)
    : detail::Extendable<conditions...>(extensions)
  {
  }

  typedef detail::boost_set_as_tparams_t<ObserverList, observers>
      observer_list_type;
  using detail::Extendable<conditions...>::get;

  template <typename input, typename result_t>
  bool
  operator()(const result_t& r, input& cache) const
  {
    // clang-format off
    static_assert(detail::all_of_v<detail::abort_condition_signature_check_v<conditions, input>...>,
                  "not all abort conditions support the specified input");
    // clang-format on

    // the different
    double conditions_stepMax = cache.step_size;
    bool   abort              = detail::abort_list_impl<conditions...>::check(
        tuple(), r, cache, conditions_stepMax);
    cache.step_size = std::min(cache.step_size, conditions_stepMax);

    return abort;
  }
};

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_HPP
