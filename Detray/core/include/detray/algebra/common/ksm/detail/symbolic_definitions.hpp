// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// This file defines a simple struct that models symbolic operations.
namespace detray::ksm::detail {
/// @brief Symbolic operators: addition, subtraction, multiplication, and
/// division.
enum class SymbolicOp { ADD, SUB, MUL, DIV };

/// @brief Symbolic operator that models A OP B.
template <SymbolicOp Op, typename A, typename B>
struct symbolic {
  static constexpr SymbolicOp operation = Op;
  using operandA = A;
  using operandB = B;
};
}  // namespace detray::ksm::detail
