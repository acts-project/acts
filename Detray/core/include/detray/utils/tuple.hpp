// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

/// Default tuple type
///
/// Serving as the final node in the recursive implementation of this tuple
/// type.
///
template <typename... Ts>
struct tuple {
  // As long as we did everything correctly, this should only get instantiated
  // with an empty parameter list, for the implementation to work correctly.
  static_assert(sizeof...(Ts) == 0, "There's a coding error in detray::tuple!");
};

/// Recursive tuple node holding one value and the remaining tuple tail.
template <typename T, typename... Ts>
struct tuple<T, Ts...> {
  // Deliberately not value-initializing: that would require every element to
  // be default constructible, which is not a requirement of this tuple type.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
  /// Default construct the recursive tuple node when all elements allow it.
  constexpr tuple() = default;

  /// Copy construct this tuple when every element is copy constructible.
  /// @param o tuple to copy
  constexpr tuple(const tuple &o)
    requires(std::is_copy_constructible_v<T> &&
             (std::is_copy_constructible_v<Ts> && ...))
  = default;

  /// Construct from a tuple with element types convertible to this tuple.
  /// @tparam U type of the first source element
  /// @tparam Us types of the remaining source elements
  /// @param o tuple to convert from
  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...>
  DETRAY_HOST_DEVICE explicit constexpr tuple(const tuple<U, Us...> &o)
      : v(o.v), r(o.r) {}

  /// Move construct this tuple when every element is movable.
  /// @param o tuple to move from
  constexpr tuple(tuple &&o) noexcept
    requires(std::is_move_constructible_v<T> &&
             (std::is_move_constructible_v<Ts> && ...))
  = default;

  /// Move construct from a tuple with compatible element types.
  /// @tparam U type of the first source element
  /// @tparam Us types of the remaining source elements
  /// @param o tuple to move from
  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...>
  DETRAY_HOST_DEVICE explicit constexpr tuple(tuple<U, Us...> &&o)
      : v(std::move(o.v)), r(std::move(o.r)) {}

  /// Construct the tuple elements from individual values.
  /// @tparam U type of the first value
  /// @tparam Us types of the remaining values
  /// @param _v value stored in this tuple node
  /// @param _r values forwarded to the tuple tail
  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...> &&
                 (!(std::is_same_v<tuple, U> ||
                    (std::is_same_v<tuple, Us> || ...)))
  DETRAY_HOST_DEVICE explicit constexpr tuple(U &&_v, Us &&..._r)
      : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

  constexpr ~tuple() noexcept = default;

  /// Copy assign this tuple when all elements are copy assignable.
  /// @param other tuple to copy
  /// @return reference to this tuple
  constexpr tuple &operator=(const tuple &other)
    requires(std::is_copy_assignable_v<T> &&
             (std::is_copy_assignable_v<Ts> && ...))
  = default;

  /// Copy assign from a tuple with compatible element types.
  /// @tparam U type of the first source element
  /// @tparam Us types of the remaining source elements
  /// @param other tuple to copy
  /// @return reference to this tuple
  template <typename U, typename... Us>
  DETRAY_HOST_DEVICE constexpr tuple &operator=(const tuple<U, Us...> &other)
    requires(std::is_assignable_v<T &, const U &> &&
             (std::is_assignable_v<Ts &, const Us &> && ...))
  {
    v = other.v;
    r = other.r;
    return *this;
  }

  /// Move assign this tuple when all elements are movable.
  /// @param other tuple to move from
  /// @return reference to this tuple
  constexpr tuple &operator=(tuple &&other) noexcept
    requires(std::is_move_assignable_v<T> &&
             (std::is_move_assignable_v<Ts> && ...))
  = default;

  /// Move assign from a tuple with compatible element types.
  /// @tparam U type of the first source element
  /// @tparam Us types of the remaining source elements
  /// @param other tuple to move from
  /// @return reference to this tuple
  template <typename U, typename... Us>
  DETRAY_HOST_DEVICE constexpr tuple &operator=(tuple<U, Us...> &&other)
    requires(std::is_assignable_v<T &, U> &&
             (std::is_assignable_v<Ts &, Us> && ...))
  {
    v = std::move(other.v);
    r = std::move(other.r);
    return *this;
  }

  T v;  ///< Value stored in this tuple node.
  // Value-initialization would require every tail element to be default
  // constructible, which is not a requirement of this tuple type.
  /// Remaining tuple elements.
  tuple<Ts...> r;  // NOLINT(cppcoreguidelines-pro-type-member-init)
};

/// A two-element tuple alias.
template <typename T1, typename T2>
using pair = tuple<T1, T2>;

/// Get a tuple element by index from a const tuple.
/// @tparam I Zero-based element index.
/// @param t Tuple to access.
/// @returns A const reference to the selected element.
template <std::size_t I, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const detray::tuple<Ts...> &t) noexcept {
  static_assert(I < sizeof...(Ts),
                "Attempt to access index greater than tuple size.");

  if constexpr (I == 0) {
    return t.v;
  } else {
    return ::detray::get<I - 1>(t.r);
  }
}

/// Get a tuple element by index from a mutable tuple.
/// @tparam I Zero-based element index.
/// @param t Tuple to access.
/// @returns A mutable reference to the selected element.
template <std::size_t I, typename... Ts>
DETRAY_HOST_DEVICE auto &get(detray::tuple<Ts...> &t) noexcept {
  static_assert(I < sizeof...(Ts),
                "Attempt to access index greater than tuple size.");

  if constexpr (I == 0) {
    return t.v;
  } else {
    return ::detray::get<I - 1>(t.r);
  }
}

/// Get a tuple element by type from a const tuple.
/// @tparam U Type of the requested element.
/// @param t Tuple containing the element.
/// @returns A const reference to the selected element.
template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const detray::tuple<T, Ts...> &t) noexcept {
  if constexpr (std::is_same_v<U, T>) {
    return t.v;
  } else if constexpr (sizeof...(Ts) > 0) {
    static_assert((std::is_same_v<U, Ts> || ...), "Type not found in tuple.");
    return ::detray::get<U, Ts...>(t.r);
  }
}

/// Get a tuple element by type from a mutable tuple.
/// @tparam U Type of the requested element.
/// @param t Tuple containing the element.
/// @returns A mutable reference to the selected element.
template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE auto &get(detray::tuple<T, Ts...> &t) noexcept {
  if constexpr (std::is_same_v<U, T>) {
    return t.v;
  } else if constexpr (sizeof...(Ts) > 0) {
    static_assert((std::is_same_v<U, Ts> || ...), "Type not found in tuple.");
    return ::detray::get<U, Ts...>(t.r);
  }
}

/// Create a tuple of references to the supplied objects.
/// @param args Objects whose references are stored.
/// @returns A tuple containing references to the arguments.
template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<Ts &...> tie(Ts &...args) {
  return detray::tuple<Ts &...>(args...);
}

/// Construct a tuple by decaying and forwarding its arguments.
/// @param vs Values to store in the tuple.
/// @returns A tuple containing the decayed argument types.
template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<std::decay_t<Ts>...> make_tuple(
    Ts &&...vs) {
  return detray::tuple<std::decay_t<Ts>...>(std::forward<Ts>(vs)...);
}

/// Construct a pair by decaying and forwarding its two arguments.
/// @param v1 First value to store.
/// @param v2 Second value to store.
/// @returns A pair containing the decayed argument types.
template <typename T1, typename T2>
DETRAY_HOST_DEVICE constexpr detray::pair<std::decay_t<T1>, std::decay_t<T2>>
make_pair(T1 &&v1, T2 &&v2) {
  return detray::pair<std::decay_t<T1>, std::decay_t<T2>>(std::forward<T1>(v1),
                                                          std::forward<T2>(v2));
}
}  // namespace detray
