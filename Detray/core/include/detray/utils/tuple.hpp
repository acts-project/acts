// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

template <typename... Ts>
struct tuple {};

template <typename T, typename... Ts>
struct tuple<T, Ts...> {
  constexpr tuple() = default;

  constexpr tuple(const tuple &o)
    requires(std::is_copy_constructible_v<T> &&
             (std::is_copy_constructible_v<Ts> && ...))
  = default;

  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...>
  DETRAY_HOST_DEVICE explicit constexpr tuple(const tuple<U, Us...> &o)
      : v(o.v), r(o.r) {}

  constexpr tuple(tuple &&o) noexcept
    requires(std::is_move_constructible_v<T> &&
             (std::is_move_constructible_v<Ts> && ...))
  = default;

  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...>
  DETRAY_HOST_DEVICE explicit constexpr tuple(tuple<U, Us...> &&o)
      : v(std::move(o.v)), r(std::move(o.r)) {}

  template <typename U, typename... Us>
    requires std::is_constructible_v<T, U &&> &&
                 std::is_constructible_v<tuple<Ts...>, Us &&...> &&
                 (!(std::is_same_v<tuple, U> ||
                    (std::is_same_v<tuple, Us> || ...)))
  DETRAY_HOST_DEVICE explicit constexpr tuple(U &&_v, Us &&..._r)
      : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

  constexpr ~tuple() noexcept = default;

  constexpr tuple &operator=(const tuple &other)
    requires(std::is_copy_assignable_v<T> &&
             (std::is_copy_assignable_v<Ts> && ...))
  = default;

  template <typename U, typename... Us>
  DETRAY_HOST_DEVICE constexpr tuple &operator=(const tuple<U, Us...> &other)
    requires(std::is_assignable_v<T &, const U &> &&
             (std::is_assignable_v<Ts &, const Us &> && ...))
  {
    v = other.v;
    r = other.r;
    return *this;
  }

  constexpr tuple &operator=(tuple &&other) noexcept
    requires(std::is_move_assignable_v<T> &&
             (std::is_move_assignable_v<Ts> && ...))
  = default;

  template <typename U, typename... Us>
  DETRAY_HOST_DEVICE constexpr tuple &operator=(tuple<U, Us...> &&other)
    requires(std::is_assignable_v<T &, U> &&
             (std::is_assignable_v<Ts &, Us> && ...))
  {
    v = std::move(other.v);
    r = std::move(other.r);
    return *this;
  }

  T v;
  tuple<Ts...> r;
};

template <typename T1, typename T2>
using pair = tuple<T1, T2>;

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

template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE const auto &get(const detray::tuple<T, Ts...> &t) noexcept {
  if constexpr (std::is_same_v<U, T>) {
    return t.v;
  } else if constexpr (sizeof...(Ts) > 0) {
    static_assert((std::is_same_v<U, Ts> || ...), "Type not found in tuple.");
    return ::detray::get<U, Ts...>(t.r);
  }
}

template <typename U, typename T, typename... Ts>
DETRAY_HOST_DEVICE auto &get(detray::tuple<T, Ts...> &t) noexcept {
  if constexpr (std::is_same_v<U, T>) {
    return t.v;
  } else if constexpr (sizeof...(Ts) > 0) {
    static_assert((std::is_same_v<U, Ts> || ...), "Type not found in tuple.");
    return ::detray::get<U, Ts...>(t.r);
  }
}

template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<Ts &...> tie(Ts &...args) {
  return detray::tuple<Ts &...>(args...);
}

template <typename... Ts>
DETRAY_HOST_DEVICE constexpr detray::tuple<std::decay_t<Ts>...> make_tuple(
    Ts &&...vs) {
  return detray::tuple<std::decay_t<Ts>...>(std::forward<Ts>(vs)...);
}

template <typename T1, typename T2>
DETRAY_HOST_DEVICE constexpr detray::pair<std::decay_t<T1>, std::decay_t<T2>>
make_pair(T1 &&v1, T2 &&v2) {
  return detray::pair<std::decay_t<T1>, std::decay_t<T2>>(std::forward<T1>(v1),
                                                          std::forward<T2>(v2));
}
}  // namespace detray
