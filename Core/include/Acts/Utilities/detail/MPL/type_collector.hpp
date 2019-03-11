// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/hana/filter.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/type.hpp>

namespace hana = boost::hana;

namespace Acts {

namespace detail {

  struct result_type_extractor
  {
    static constexpr auto predicate = hana::is_valid(
        [](auto t) -> hana::type<typename decltype(t)::type::result_type>{});
    template <typename T>
    using extractor_impl            = typename T::result_type;
    static constexpr auto extractor = hana::template_<extractor_impl>;
  };

  struct action_type_extractor
  {
    static constexpr auto predicate = hana::is_valid(
        [](auto t) -> hana::type<typename decltype(t)::type::action_type>{});
    template <typename T>
    using extractor_impl            = typename T::action_type;
    static constexpr auto extractor = hana::template_<extractor_impl>;
  };

  constexpr auto type_collector = [](auto t_, auto predicate, auto extractor) {
    constexpr auto have_result
        = hana::filter(t_, [&](auto t) { return predicate(t); });
    constexpr auto result_types
        = hana::to_set(hana::transform(have_result, extractor));
    return result_types;
  };

  template <typename helper, typename... traits>
  constexpr auto type_collector_t = type_collector(hana::tuple_t<traits...>,
                                                   helper::predicate,
                                                   helper::extractor);

  template <typename T>
  constexpr bool has_result_type_v
      = decltype(result_type_extractor::predicate(hana::type_c<T>))::value;

  template <typename T>
  using result_type_t = typename result_type_extractor::extractor_impl<T>;

  template <typename T>
  constexpr bool has_action_type_v
      = decltype(action_type_extractor::predicate(hana::type_c<T>))::value;

  template <typename T>
  using action_type_t = typename action_type_extractor::extractor_impl<T>;
}  // namespace detail

}  // namespace Acts
