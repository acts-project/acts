// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/action_list_implementation.hpp"
#include "Acts/Propagator/detail/action_signature_check.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

/// @brief ActionList implementation to be used with the propagator
///
/// This is the ActionList struct that is used in the propagator
/// to define a list of different actors_t that are eacch
/// executed during the stepping procedure
template <typename... actors_t>
struct ActionList : public detail::Extendable<actors_t...>
{
private:
  static_assert(not detail::has_duplicates_v<actors_t...>,
                "same action type specified several times");

  // clang-format off
  using results = detail::type_collector_t<detail::result_type_extractor, actors_t...>;
  // clang-format on

  using detail::Extendable<actors_t...>::tuple;

public:
  template <template <typename...> class R>
  using result_type = detail::boost_set_as_tparams_t<R, results>;

  using detail::Extendable<actors_t...>::get;

  /// Default constructor
  ActionList() = default;

  /// Default copy constructor
  ///
  /// @param actors The source action list
  ActionList(const ActionList<actors_t...>& actors) = default;

  /// Default move constructor
  ///
  /// @param actors The source action list
  ActionList(ActionList<actors_t...>&& actors) = default;

  /// Default move assignment operator
  ///
  /// @param actors The source action list
  ActionList<actors_t...>&
  operator=(const ActionList<actors_t...>& actors)
      = default;

  /// Default move assignment operator
  ///
  /// @param actors The source action list
  ActionList<actors_t...>&
  operator=(ActionList<actors_t...>&& actors)
      = default;

  /// Call operator that is that broadcasts the call to the tuple()
  /// members of the list
  ///
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper used for the propagation
  /// @tparam result_t is the result type from actions
  ///
  /// @param [in,out] prop is the propagator state object
  /// @param [in] stepper The stepper in use
  /// @param [in,out] result is the result object from actions
  ///
  /// @return bool type indiciating if the step size can be released
  template <typename propagator_state_t, typename stepper_t, typename result_t>
  void
  operator()(propagator_state_t& state,
             const stepper_t&    stepper,
             result_t&           result) const
  {
    // clang-format off
    static_assert(detail::all_of_v<detail::action_signature_check_v<actors_t, 
                                      propagator_state_t, stepper_t>...>,
                  "not all actors support the method signature");
    // clang-format on

    using impl = detail::action_list_impl<actors_t...>;
    impl::action(tuple(), state, stepper, result);
  }
};

}  // namespace Acts