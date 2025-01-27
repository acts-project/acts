// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <type_traits>

namespace Acts {

/// @brief type list implementation
/// @see https://www.codingwiththomas.com/blog/getting-started-with-typelists
template <typename... Ts>
struct TypeList {};

namespace Types {

/// Number of types in the list
/// @{
template <typename = void>
struct getSize {};

template <typename... Ts>
struct getSize<TypeList<Ts...>>
    : std::integral_constant<std::size_t, sizeof...(Ts)> {};

template <typename L>
constexpr inline std::size_t size{getSize<L>()};
/// @}

/// Access the first type
/// @{
template <typename = void>
struct getFront {};

template <typename T, typename... Ts>
struct getFront<TypeList<T, Ts...>> {
  using type = T;
};

template <typename L>
using front = typename getFront<L>::type;
/// @}

/// Access the last type
/// @{
template <typename = void>
struct getBack {};

template <typename T, typename... Ts>
struct getBack<TypeList<T, Ts...>> {
  using type = std::conditional_t<sizeof...(Ts) == 0, T,
                                  typename getBack<TypeList<Ts...>>::type>;
};
// Base case
template <>
struct getBack<TypeList<>> {
  using type = void;
};
template <typename L>
using back = typename getBack<L>::type;
/// @}

/// Append a type
/// @{
template <typename N, typename = void>
struct doPushBack {};

template <typename N, typename... Ts>
struct doPushBack<N, TypeList<Ts...>> {
  using type = TypeList<Ts..., N>;
};
template <typename L, typename N>
using push_back = typename doPushBack<N, L>::type;
/// @}

/// Prepend a type
/// @{
template <typename N, typename = void>
struct doPushFront {};

template <typename N, typename... Ts>
struct doPushFront<N, TypeList<Ts...>> {
  using type = TypeList<N, Ts...>;
};

template <typename L, typename N>
using push_front = typename doPushFront<N, L>::type;
/// @}

}  // namespace Types
}  // namespace Acts
