// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <cstddef>

// This file describes the values that a KSM matrix can take.
namespace detray::ksm {
/// @brief A constant integral value.
///
/// Note we don't support floating point constants right now, as static
/// floating points can become a mess in terms of correctness. Two constants
/// that are equal but not bit-identical would be distinct types, so
/// canonicalization could never unify them; and eliding x * 0 is wrong when x
/// is NaN or infinite.
template <int I>
struct integral_value {
  static constexpr int value = I;
};

/// @brief An anonymous variable.
struct variable {};

/// @brief A non-anonymous variable.
template <std::size_t I>
struct index_variable {
  static constexpr std::size_t index = I;
};

/// @brief A constant value zero.
using zero = integral_value<0>;

/// @brief A constant value one.
using one = integral_value<1>;
}  // namespace detray::ksm
