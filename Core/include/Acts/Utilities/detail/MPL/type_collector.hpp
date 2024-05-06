// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/hana/filter.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/type.hpp>

namespace Acts::detail {
namespace hana = boost::hana;

/**
 * Struct which extracts the result type from an actor.
 * This is used as an argument to `type_collector`.
 */
struct result_type_extractor {
  // Checks whether the type even has a result type
  static constexpr auto predicate = hana::is_valid(
      [](auto t) -> hana::type<typename decltype(t)::type::result_type> {});
  // meta function to extract the result type
  template <typename T>
  using extractor_impl = typename T::result_type;
  // embed meta function in hana
  static constexpr auto extractor = hana::template_<extractor_impl>;
};

/**
 * Struct which extracts the action type from an aborter.
 * This is used as an argument to `type_collector`.
 */
struct action_type_extractor {
  // Checks if aborter even has action type
  static constexpr auto predicate = hana::is_valid(
      [](auto t) -> hana::type<typename decltype(t)::type::action_type> {});
  // meta function to extract the action type
  template <typename T>
  using extractor_impl = typename T::action_type;
  // converted to hana
  static constexpr auto extractor = hana::template_<extractor_impl>;
};

/**
 * The main type collector. This loops over the given tuple of actions or
 * aborters,
 * filters by predicate and uses extractor to construct a resulting output
 * set.
 */
constexpr auto type_collector = [](auto t_, auto predicate, auto extractor) {
  // filtered list using predicate
  constexpr auto have_result =
      hana::filter(t_, [&](auto t) { return predicate(t); });
  // convert to set to remove duplicates, and transform to unpacked type
  // using extractor.
  constexpr auto result_types =
      hana::to_set(hana::transform(have_result, extractor));
  return result_types;
};

/**
 * Helper around type_collector which constructrs a hana tuple input from
 * variadic
 * template args, and pre-unpacks the predicate and extractor from the helper
 * type
 * @tparam helper Either result_type_extractor or action_type_extractor
 * @tparam items The items to filter / collect from.
 */
template <typename helper, typename... items>
constexpr auto type_collector_t = type_collector(hana::tuple_t<items...>,
                                                 helper::predicate,
                                                 helper::extractor);

/**
 * Meta function which returns a compile time bool
 * @tparam T the type to check
 */
template <typename T>
constexpr bool has_result_type_v =
    decltype(result_type_extractor::predicate(hana::type_c<T>))::value;

/**
 * Meta function which gets the result type from an action
 * @tparam T The type to extract from
 */
template <typename T>
using result_type_t = typename result_type_extractor::extractor_impl<T>;

/**
 * Meta function which returns a compile time bool
 * @tparam T the type to check
 */
template <typename T>
constexpr bool has_action_type_v =
    decltype(action_type_extractor::predicate(hana::type_c<T>))::value;

/**
 * Meta function which gets the action for an aborter
 * @tparam T The type to extract from
 */
template <typename T>
using action_type_t = typename action_type_extractor::extractor_impl<T>;
}  // namespace Acts::detail
