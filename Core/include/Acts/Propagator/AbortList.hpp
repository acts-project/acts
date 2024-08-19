// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/abort_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <boost/hana/type.hpp>
#include <boost/hana/unpack.hpp>

namespace hana = boost::hana;

namespace Acts {

/// @brief AbortList object to be used in the propagation
///
/// The abort list is a list of structs or classes that
/// is called at each propagation step and can trigger the abort
/// of the current propagation.
///
/// It can (optionally) depend on a result of an Actor from
/// the actor list.
template <typename... aborters_t>
struct AbortList : public detail::Extendable<aborters_t...> {
 private:
  static_assert(!detail::has_duplicates_v<aborters_t...>,
                "same aborter type specified several times");

  using detail::Extendable<aborters_t...>::tuple;

 public:
  /// @cond
  // This uses the type collector
  using result_type = typename decltype(hana::unpack(
      detail::type_collector_t<detail::action_type_extractor, aborters_t...>,
      hana::template_<AbortList>))::type;
  /// @endcond

  using detail::Extendable<aborters_t...>::get;

  /// Default constructor
  AbortList() = default;

  /// Default copy constructor
  ///
  /// @param aborters The source action list
  AbortList(const AbortList<aborters_t...>& aborters) = default;

  /// Default move constructor
  ///
  /// @param aborters The source action list
  AbortList(AbortList<aborters_t...>&& aborters) = default;

  /// Default move assignment operator
  ///
  /// @param aborters The source action list
  AbortList<aborters_t...>& operator=(
      const AbortList<aborters_t...>& aborters) = default;

  /// Default move assignment operator
  ///
  /// @param aborters The source action list
  AbortList<aborters_t...>& operator=(AbortList<aborters_t...>&& aborters) =
      default;

  /// Constructor from tuple
  ///
  /// @param aborters Source extensions tuple
  AbortList(const std::tuple<aborters_t...>& aborters)
      : detail::Extendable<aborters_t...>(aborters) {}

  /// Constructor from tuple move
  ///
  /// @param aborters Source extensions tuple
  AbortList(std::tuple<aborters_t...>&& aborters)
      : detail::Extendable<aborters_t...>(std::move(aborters)) {}

  /// Append new entries and return a new condition
  template <typename... appendices_t>
  AbortList<aborters_t..., appendices_t...> append(appendices_t... aps) const {
    auto catTuple =
        std::tuple_cat(tuple(), std::tuple<appendices_t...>(aps...));
    return AbortList<aborters_t..., appendices_t...>(std::move(catTuple));
  }

  /// This is the call signature for the abort list, it broadcasts the call
  /// to the tuple() members of the list
  ///
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state is the state object from the propagator
  /// @param [in] stepper Stepper used for the propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param [in] args are the arguments to be passed to the aborters
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  bool operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, Args&&... args) const {
    using impl = detail::abort_list_impl<aborters_t...>;
    return impl::check(tuple(), state, stepper, navigator,
                       std::forward<Args>(args)...);
  }
};

}  // namespace Acts
