// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts {
/// @cond detail
namespace detail {
/**
 * @brief get position of integral constant in template parameter pack
 *
 * @tparam T integral type of the values to be investigated
 * @tparam target target value whose position in the template parameter pack
 * should be determined
 * @tparam values template parameter pack containing the list of values
 *
 * @return `get_position<T,target,values...>::value` yields the position of
 * `target` inside `values`.
 *         If `target` is not in the list of `values`, a compile-time error is
 * generated.
 */
template <typename T, T target, T... values>
struct get_position;

/// @cond
template <typename T, T target, T... others>
struct get_position<T, target, target, others...> {
  enum { value = 0 };
};

template <typename T, T target, T next, T... others>
struct get_position<T, target, next, others...> {
  enum { value = get_position<T, target, others...>::value + 1 };
};
/// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts
