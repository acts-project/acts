// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ACTION_LIST_HPP
#define ACTS_ACTION_LIST_HPP

#include "ACTS/Propagator/detail/Extendable.hpp"
#include "ACTS/Propagator/detail/action_list_implementation.hpp"
#include "ACTS/Propagator/detail/action_signature_check.hpp"
#include "ACTS/Utilities/detail/MPL/all_of.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

/// This is the ActionList struct that is used in the propagator
/// to define a list of different actions that are eacch
/// executed during the stepping procedure
template <typename... actions>
struct ActionList : private detail::Extendable<actions...>
{
private:
  static_assert(not detail::has_duplicates_v<actions...>,
                "same action type specified several times");

  // clang-format off
  typedef detail::type_collector_t<detail::result_type_extractor, actions...> results;
  // clang-format on

  using detail::Extendable<actions...>::tuple;

public:
  template <template <typename...> class R>
  using result_type = detail::boost_set_as_tparams_t<R, results>;

  using detail::Extendable<actions...>::get;

  /// Call operator that is that broadcasts the call to the tuple()
  /// members of the list
  ///
  /// @tparam cache_t is the cache type from the propagator
  /// @tparam result_t is the result type from actions
  ///
  /// @param cache[in,out] is the cache object from the propagator
  /// @param result[in,out] is the result object from actions
  ///
  /// @return bool type indiciating if the step size can be released
  template <typename cache_t, typename result_t>
  void
  operator()(cache_t& cache, result_t& result) const
  {
    // clang-format off
    static_assert(detail::all_of_v<detail::action_signature_check_v<actions, cache_t>...>,
                  "not all actions support the specified cache_t");
    // clang-format on

    typedef detail::action_list_impl<actions...> impl;
    impl::action(tuple(), cache, result);
  }
};

}  // namespace Acts

#endif  // ACTS_ACTION_LIST_HPP
