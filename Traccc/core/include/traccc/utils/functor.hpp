// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::details::functor {
/**
 * @brief The identity functor, such that `identity<T>` is equivalent to `T`.
 */
template <typename T>
using identity = T;

/**
 * @brief The reference functor, such that `reference<T>` is equivalent to `T&`.
 */
template <typename T>
using reference = T&;

/**
 * @brief The const reference functor, such that `const_reference<T>` is
 * equivalent to `const T&`.
 */
template <typename T>
using const_reference = const T&;

/**
 * @brief A natural transformation between two functors.
 *
 * Given two functors F1 and F2 as well as some types Ts... such that
 * `F2<Ts...>` is a type, this produces `F1<Ts...>`.
 */
template <template <typename...> typename F, typename T>
struct reapply {};

template <template <typename...> typename F1,
          template <typename...> typename F2, typename... Ts>
struct reapply<F1, F2<Ts...>> {
  using type = F1<Ts...>;
};
}  // namespace traccc::details::functor
