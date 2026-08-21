// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"
#include "detray/algebra/common/ksm/detail/symbolic_definitions.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <concepts>

// This file describes how symbolic operations are built over values, and how
// they are resolved back down to a value.
namespace detray::ksm::detail {
/// @brief Compute the negation of a value.
template <typename V>
constexpr auto negate_value() {
  if constexpr (concepts::is_variable<V>) {
    // The negation of a variable, is simply another variable.
    static_assert(concepts::is_index_variable<V>);
    return V{};
  } else {
    // The negation of a constant is also just a constant with the sign
    // flipped.
    return integral_value<-V::value>{};
  }
}

/// @brief Resolve a symbolic operation.
///
/// @tparam Symbolic The operation to resolve.
/// @tparam FromSameMatrix Whether all operands in the operation are from the
/// same matrix.
/// @tparam CollapseToVariable Whether to collapse values to variables, or to
/// keep them in simplified symbolic form.
template <typename Symbolic, bool FromSameMatrix = false,
          bool CollapseToVariable = true>
constexpr auto resolve_symbolic() {
  // Check whether we are performing a division by zero statically.
  static_assert(!(Symbolic::operation == SymbolicOp::DIV &&
                  std::same_as<typename Symbolic::operandB, zero>));

  // This code now branches (statically) across all possible operation types.
  if constexpr (Symbolic::operation == SymbolicOp::ADD) {
    // For addition, we want to model the following function, where o is an
    // arbitrary operand; if collapsing to a variable is turned on, the result
    // at the end between parentheses is returned, if different from the result
    // otherwise:
    //
    // (+) (Variable i) o = Symbolic + (Variable i) o (= AnonVariable)
    // (+) o (Variable i) = Symbolic + o (Variable i) (= AnonVariable)
    // (+) (Constant i) (Constant j) = Constant (i + j)
    // (+) oi oj = Symbolic + oi oj (= AnonVariable)
    //
    // Note that the last line is unnecessarily generic; if oi or oj are
    // symbolic operations, we can simplify further. But for now we just
    // collapse them.
    if constexpr (concepts::is_variable<typename Symbolic::operandA> ||
                  concepts::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else if constexpr (concepts::is_integral_value<
                             typename Symbolic::operandA> &&
                         concepts::is_integral_value<
                             typename Symbolic::operandB>) {
      return integral_value<Symbolic::operandA::value +
                            Symbolic::operandB::value>{};
    } else {
      // Neither operand is a variable, and they are not both constants, so at
      // least one of them has to be a symbolic operation.
      static_assert(concepts::is_symbolic<typename Symbolic::operandA> ||
                    concepts::is_symbolic<typename Symbolic::operandB>);
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::SUB) {
    // For subtraction, we model an almost identical function as addition:
    //
    // (-) (Variable i) o = Symbolic - (Variable i) o (= AnonVariable)
    // (-) o (Variable i) = Symbolic - o (Variable i) (= AnonVariable)
    // (-) (Constant i) (Constant j) = Constant (i - j)
    // (-) oi oj = Symbolic - oi oj (= AnonVariable)
    //
    // With the critical difference that, if the objects are from the same
    // matrix, we can return zero if subtracting a variable from itself.
    if constexpr (FromSameMatrix && std::same_as<typename Symbolic::operandA,
                                                 typename Symbolic::operandB>) {
      return zero{};
    } else if constexpr (concepts::is_variable<typename Symbolic::operandA> ||
                         concepts::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else if constexpr (concepts::is_integral_value<
                             typename Symbolic::operandA> &&
                         concepts::is_integral_value<
                             typename Symbolic::operandB>) {
      return integral_value<Symbolic::operandA::value -
                            Symbolic::operandB::value>{};
    } else {
      // As with addition: at least one operand must be a symbolic operation.
      static_assert(concepts::is_symbolic<typename Symbolic::operandA> ||
                    concepts::is_symbolic<typename Symbolic::operandB>);
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::MUL) {
    // For multiplication we model the following:
    //
    // (*) o (Constant 0) = Constant 0
    // (*) (Constant 0) o = Constant 0
    // (*) (Variable i) o = Symbolic * (Variable i) o (= AnonVariable)
    // (*) o (Variable i) = Symbolic * o (Variable i) (= AnonVariable)
    // (*) (Constant i) (Constant j) = Constant (i * j)
    // (*) oi oj = Symbolic * oi oj (= AnonVariable)
    if constexpr (std::same_as<typename Symbolic::operandA, zero> ||
                  std::same_as<typename Symbolic::operandB, zero>) {
      return zero{};
    } else if constexpr (concepts::is_variable<typename Symbolic::operandA> ||
                         concepts::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else if constexpr (concepts::is_integral_value<
                             typename Symbolic::operandA> &&
                         concepts::is_integral_value<
                             typename Symbolic::operandB>) {
      return integral_value<Symbolic::operandA::value *
                            Symbolic::operandB::value>{};
    } else {
      // As with addition: at least one operand must be a symbolic operation.
      static_assert(concepts::is_symbolic<typename Symbolic::operandA> ||
                    concepts::is_symbolic<typename Symbolic::operandB>);
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::DIV) {
    // For division we model the most complex function, but we simplify it
    // quite a bit:
    //
    // (/) o (Constant 0) = Error
    // (/) o (Constant 1) = o
    // (/) (Constant 0) o = Constant 0
    // (/) oi oj = Symbolic / oi oj (= AnonVariable)
    if constexpr (std::same_as<typename Symbolic::operandB, one>) {
      return typename Symbolic::operandA{};
    } else if constexpr (std::same_as<typename Symbolic::operandA, zero>) {
      return zero{};
    } else {
      // Could handle division of constants here, but has only marginal returns.
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  }
}

/// @brief Resolver for any value, including ones that are not symbolic
/// operations.
///
/// @{
template <typename V, bool = false, bool = true>
struct resolve_value {};

// For symbolic operations, call into resolve_symbolic.
template <typename V, bool FromSameMatrix, bool CollapseToVariable>
  requires(concepts::is_symbolic<V>)
struct resolve_value<V, FromSameMatrix, CollapseToVariable> {
  using type =
      decltype(resolve_symbolic<V, FromSameMatrix, CollapseToVariable>());
};

// For non-operations, simply return the input!
template <typename V, bool FromSameMatrix, bool CollapseToVariable>
  requires(!concepts::is_symbolic<V>)
struct resolve_value<V, FromSameMatrix, CollapseToVariable> {
  using type = V;
};
/// @}

/// @brief Helper functions for common operations
///
/// @{
template <typename A, typename B>
  requires(concepts::is_value<A> && concepts::is_value<B>)
using value_addition = symbolic<SymbolicOp::ADD, A, B>;

template <typename A, typename B>
  requires(concepts::is_value<A> && concepts::is_value<B>)
using value_subtraction = symbolic<SymbolicOp::SUB, A, B>;

template <typename A, typename B>
  requires(concepts::is_value<A> && concepts::is_value<B>)
using value_multiplication = symbolic<SymbolicOp::MUL, A, B>;

template <typename A, typename B>
  requires(concepts::is_value<A> && concepts::is_value<B>)
using value_division = symbolic<SymbolicOp::DIV, A, B>;

template <typename V>
  requires(concepts::is_value<V>)
using value_negation = decltype(negate_value<V>());
/// @}
}  // namespace detray::ksm::detail
