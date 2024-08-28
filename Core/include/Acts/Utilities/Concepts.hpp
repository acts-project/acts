// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <concepts>
#include <type_traits>

namespace Acts::Concepts {
/// @brief Concept that is true if T is the same as any of Ts.
template <typename T, typename... Ts>
concept same_as_any_of = (std::same_as<T, Ts> || ...);

/// @brief Concept that is equivalent to `is_nothrow_move_constructible`.
/// @todo Convert this to a "real" concept.
template <typename T>
concept nothrow_move_constructible =
    std::is_nothrow_move_constructible<T>::value;

/// @brief Concept that is true if T is an arithmetic type.
template <typename T>
concept arithmetic = std::integral<T> || std::floating_point<T>;

/// @brief Concept that is satisfied iff both of its arguments decay to the
/// same type.
template <typename T1, typename T2>
concept decayed_same_as = std::same_as<std::decay_t<T1>, std::decay_t<T2> >;
}  // namespace Acts::Concepts
