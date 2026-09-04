// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <concepts>
#include <type_traits>

namespace detray::concepts {

template <typename T, typename... U>
concept any_of = std::disjunction_v<std::is_same<T, U>...>;

/// Same, except for cv qualifiers and reference
template <typename T, typename U>
concept same_as_cvref =
    std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;

/// Concept for detecting when a type is a non-const version of another
template <typename T, typename U>
concept same_as_no_const = std::same_as<std::remove_cv_t<T>, U>;

/// Type identifier concept
template <typename T>
concept type_id = std::is_enum_v<T>;

/// Concept that checks if a type models an interval of some value that can
/// be obtained with 'get'.
template <typename I>
concept interval = requires(I i) {
  requires(!std::is_fundamental_v<std::remove_cvref_t<I>>);

  { detray::detail::get<0>(i) } -> concepts::arithmetic_cvref;

  { detray::detail::get<1>(i) } -> concepts::arithmetic_cvref;
};

/// Concept that checks whether a type can be incremented with arbitrary
/// integer values.
template <typename T>
concept random_access_incrementable = requires(T i, const T j, int n) {
  { i += n } -> std::same_as<T &>;
  { j + n } -> std::same_as<T>;
  { n + j } -> std::same_as<T>;
  { i -= n } -> std::same_as<T &>;
  { j - n } -> std::same_as<T>;
};
}  // namespace detray::concepts
