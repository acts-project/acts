// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <utility>

namespace Acts {

namespace detail {

/// @brief Base case: apply a single function to a value.
/// @param x The input value
/// @param f The function to apply
/// @return The result of invoking @p f with @p x
template <typename T, typename F>
decltype(auto) applyRight(T&& x, const F& f) {
  return std::invoke(f, std::forward<T>(x));
}

/// @brief Recursive case: apply the last function first, then fold right.
/// @param x The input value
/// @param f The outermost function to apply last
/// @param rest The remaining functions, applied right-to-left
/// @return The result of the right-to-left function chain
template <typename T, typename F, typename... Rest>
decltype(auto) applyRight(T&& x, const F& f, const Rest&... rest) {
  return std::invoke(f, applyRight(std::forward<T>(x), rest...));
}

}  // namespace detail

/// @brief Compose multiple callables into a single callable, applied
///        right-to-left (mathematical convention).
///
/// Given `compose(f, g, h)`, the returned callable computes `f(g(h(x)))`.
/// All callables are stored exactly once in a tuple, avoiding the nested-copy
/// overhead of a recursive lambda approach.
///
/// @tparam Fs The callable types
/// @param fs The callables to compose
/// @return A callable that applies @p fs right-to-left
template <typename... Fs>
auto compose(Fs&&... fs) {
  return [tup = std::make_tuple(std::forward<Fs>(fs)...)]<typename T>(
             T&& x) -> decltype(auto) {
    return std::apply(
        [&](const auto&... fns) -> decltype(auto) {
          return detail::applyRight(std::forward<T>(x), fns...);
        },
        tup);
  };
}

}  // namespace Acts
