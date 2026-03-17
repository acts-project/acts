// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Acts::Ranges {

namespace detail {

// Forward declaration to allow to_adaptor::operator() to return
// bound_to_adaptor
template <template <typename...> class Container, typename... BoundArgs>
struct bound_to_adaptor;

/// Adaptor for converting a range to a container
/// @tparam Container the container type to convert to
template <template <typename...> class Container>
struct to_adaptor {
  /// Convert a range to a container
  /// @tparam Range the range type to convert
  /// @tparam Args additional arguments forwarded to the container constructor
  /// @param range the range to convert
  /// @return the converted container
  template <std::ranges::input_range Range, typename... Args>
  auto operator()(Range&& range, Args&&... args) const {
    using ValueType = std::ranges::range_value_t<std::remove_cvref_t<Range>>;
    if constexpr (requires {
                    typename Container<ValueType, std::decay_t<Args>...>;
                  }) {
      return Container<ValueType, std::decay_t<Args>...>(
          std::ranges::begin(range), std::ranges::end(range),
          std::forward<Args>(args)...);
    } else {
      return Container<ValueType>(std::ranges::begin(range),
                                  std::ranges::end(range),
                                  std::forward<Args>(args)...);
    }
  }

  /// Create a bound adaptor for piping with extra constructor arguments
  /// @tparam Arg0 the first extra argument (must not itself be an input range)
  /// @tparam Rest the remaining extra argument types
  /// @return a bound_to_adaptor that can be used with operator|
  template <typename Arg0, typename... Rest>
    requires(!std::ranges::input_range<std::remove_cvref_t<Arg0>>)
  auto operator()(Arg0&& arg0, Rest&&... rest) const {
    return bound_to_adaptor<Container, std::decay_t<Arg0>,
                            std::decay_t<Rest>...>{
        std::make_tuple(std::forward<Arg0>(arg0), std::forward<Rest>(rest)...)};
  }
};

/// Adaptor for converting a range to a container with pre-bound constructor
/// arguments, enabling use with operator|
/// @tparam Container the container type to convert to
/// @tparam BoundArgs the pre-bound constructor argument types
template <template <typename...> class Container, typename... BoundArgs>
struct bound_to_adaptor {
  std::tuple<BoundArgs...> args;

  /// Convert a range to a container using the pre-bound arguments
  /// @tparam Range the range type to convert
  /// @param range the range to convert
  /// @return the converted container
  template <std::ranges::input_range Range>
  auto operator()(Range&& range) const {
    return std::apply(
        [&range](auto&&... a) {
          return to_adaptor<Container>{}(std::forward<Range>(range),
                                         std::forward<decltype(a)>(a)...);
        },
        args);
  }
};

template <typename Range, template <typename...> class Container>
auto operator|(Range&& range, to_adaptor<Container> adaptor) {
  return adaptor(std::forward<Range>(range));
}

template <typename Range, template <typename...> class Container,
          typename... BoundArgs>
auto operator|(Range&& range,
               const bound_to_adaptor<Container, BoundArgs...>& adaptor) {
  return adaptor(std::forward<Range>(range));
}

}  // namespace detail

/// Convert a range to a container
/// @tparam Container the container type to convert to
template <template <typename...> class Container>
constexpr detail::to_adaptor<Container> to{};

}  // namespace Acts::Ranges
