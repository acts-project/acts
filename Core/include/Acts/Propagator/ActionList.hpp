// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/action_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <boost/hana/type.hpp>
#include <boost/hana/unpack.hpp>

namespace hana = boost::hana;

namespace Acts {

/// @brief ActionList implementation to be used with the propagator
///
/// This is the ActionList struct that is used in the propagator
/// to define a list of different actors_t that are each
/// executed during the stepping procedure
template <typename... actors_t>
struct ActionList : public detail::Extendable<actors_t...> {
 private:
  static_assert(!detail::has_duplicates_v<actors_t...>,
                "same action type specified several times");

  using detail::Extendable<actors_t...>::tuple;

 public:
  /// @cond
  // This uses the type collector and unpacks using the `R` meta function
  template <template <typename...> class R>
  using result_type = typename decltype(hana::unpack(
      detail::type_collector_t<detail::result_type_extractor, actors_t...>,
      hana::template_<R>))::type;
  /// @endcond

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
  ActionList<actors_t...>& operator=(const ActionList<actors_t...>& actors) =
      default;

  /// Default move assignment operator
  ///
  /// @param actors The source action list
  ActionList<actors_t...>& operator=(ActionList<actors_t...>&& actors) =
      default;

  /// Call operator that is that broadcasts the call to the tuple()
  /// members of the list
  ///
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper used for the propagation
  /// @tparam navigator_t Type of the navigator used for the propagation
  ///
  /// @param [in,out] state This is the propagator state object
  /// @param [in] stepper The stepper in use
  /// @param [in] navigator The navigator in use
  /// @param [in] args The arguments to be passed to the actions
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, Args&&... args) const {
    using impl = detail::action_list_impl<actors_t...>;
    impl::action(tuple(), state, stepper, navigator,
                 std::forward<Args>(args)...);
  }
};

}  // namespace Acts
