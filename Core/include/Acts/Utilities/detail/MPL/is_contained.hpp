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
 * @brief check whether given integral constant is contained in a template
 * parameter pack
 *
 * @tparam T integral type of the values to be checked
 * @tparam target target value to be looked for
 * @tparam values template parameter pack containing the list of values
 *
 * @return `is_contained<T,target,values...>::value` is @c true if `target` is
 * among `values`, otherwise @c false
 */
template <typename T, T target, T... values>
struct is_contained;

/// @cond
template <typename T, T target, T... others>
struct is_contained<T, target, target, others...> {
  enum { value = true };
};

template <typename T, T target>
struct is_contained<T, target, target> {
  enum { value = true };
};

template <typename T, T target, T next, T... others>
struct is_contained<T, target, next, others...> {
  enum { value = is_contained<T, target, others...>::value };
};

template <typename T, T target, T next>
struct is_contained<T, target, next> {
  enum { value = false };
};
/// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts
