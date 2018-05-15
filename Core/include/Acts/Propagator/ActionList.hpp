// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/Extendable.hpp"
#include "Acts/Propagator/detail/action_list_implementation.hpp"
#include "Acts/Propagator/detail/action_signature_check.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

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
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_state_t is the state type of the stepper
  /// @tparam result_t is the result type from actions
  ///
  /// @param propState[in,out] is the propagator state object
  /// @param stepState[in,out] is the stepper state object
  /// @param result[in,out] is the result object from actions
  ///
  /// @return bool type indiciating if the step size can be released
  template <typename propagator_state_t,
            typename stepper_state_t,
            typename result_t>
  void
  operator()(propagator_state_t& propState,
             stepper_state_t&    stepState,
             result_t&           result) const
  {
    // clang-format off
    static_assert(detail::all_of_v<detail::action_signature_check_v<actions, 
                                      propagator_state_t, 
                                      stepper_state_t>...>,
                  "not all actions support the method signature");
    // clang-format on

    typedef detail::action_list_impl<actions...> impl;
    impl::action(tuple(), propState, stepState, result);
  }
};

}  // namespace Acts