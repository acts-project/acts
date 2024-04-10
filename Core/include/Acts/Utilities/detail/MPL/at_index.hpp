// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts::detail {
/**
 * @brief return integral constant at position in template parameter pack
 *
 * @tparam T integral type of the values to be investigated
 * @tparam index position in the template parameter pack
 * @tparam values template parameter pack containing the list of values
 *
 * @return `at_index<T,index,values...>::type` yields the integral constant
 *         at the given position in `values` if 0 &le; `index` <
 *         sizeof...(values). Otherwise, a compile-time error is generated.
 */
template <typename T, std::size_t index, T... values>
struct at_index;

/// @cond
template <typename T, std::size_t index, T next, T... others>
struct at_index<T, index, next, others...> {
  static constexpr T value = at_index<T, index - 1, others...>::value;
};

template <typename T, T next, T... others>
struct at_index<T, 0, next, others...> {
  static constexpr T value = next;
};
/// @endcond
}  // namespace Acts::detail
