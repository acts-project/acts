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
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/assert.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <array>
#include <concepts>
#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

/// @brief Matrices whose substructure is known at compile time.
///
/// A matrix cell is described by a type rather than a value: a structural
/// @c zero or @c one, some other compile-time @c integral_value, or a
/// @c variable that is actually stored. Cells that name the same
/// @c index_variable share a single stored scalar, which is how symmetry and
/// other repeated structure are expressed.
///
/// Operations are performed on the substructure as well as on the values, so
/// the result type of a product or a congruence knows its own structure, and
/// the emitted code skips the multiplications by zero and one, the additions
/// of constants, and the duplicate storage that a structure-blind matrix
/// would pay for.
namespace detray::ksm {

template <int I>
struct integral_value {
  static constexpr int value = I;
};
struct variable {};
template <std::size_t I>
struct index_variable {
  static constexpr std::size_t index = I;
};
using zero = integral_value<0>;
using one = integral_value<1>;

template <typename... value_ts>
struct row;

template <typename... row_ts>
struct substructure;

namespace detail {

enum class SymbolicOp { ADD, SUB, MUL, DIV };

template <SymbolicOp Op, typename A, typename B>
struct symbolic {
  static constexpr SymbolicOp operation = Op;
  using operandA = A;
  using operandB = B;
};

namespace checks {
namespace detail {
template <typename T>
struct is_integral_value_helper {
  static constexpr bool value = false;
};

template <int I>
struct is_integral_value_helper<integral_value<I>> {
  static constexpr bool value = true;
};

template <typename T>
struct is_symbolic_helper {
  static constexpr bool value = false;
};

template <SymbolicOp op, typename A, typename B>
struct is_symbolic_helper<symbolic<op, A, B>> {
  static constexpr bool value = true;
};

template <typename T>
struct is_index_variable_helper {
  static constexpr bool value = false;
};

template <std::size_t I>
struct is_index_variable_helper<index_variable<I>> {
  static constexpr bool value = true;
};
}  // namespace detail

template <typename T>
concept is_integral_value = detail::is_integral_value_helper<T>::value;

template <typename T>
concept is_symbolic = detail::is_symbolic_helper<T>::value;

template <typename T>
concept is_index_variable = detail::is_index_variable_helper<T>::value;

template <typename T>
concept is_variable = is_index_variable<T> || std::same_as<T, variable>;

template <typename T>
concept is_value = is_symbolic<T> || is_integral_value<T> || is_variable<T>;

template <typename T>
concept is_canonical_value = is_integral_value<T> || is_index_variable<T>;

namespace detail {
template <typename T>
inline constexpr bool is_row_v = false;
template <typename... Vs>
inline constexpr bool is_row_v<row<Vs...>> = (is_value<Vs> && ...);

template <typename T>
inline constexpr bool is_canonical_row_v = false;
template <typename... Vs>
inline constexpr bool is_canonical_row_v<row<Vs...>> =
    (is_canonical_value<Vs> && ...);

template <typename T>
inline constexpr bool is_canonical_substructure_v = false;
template <typename... Rs>
inline constexpr bool is_canonical_substructure_v<substructure<Rs...>> =
    (is_canonical_row_v<Rs> && ...);

template <typename Sub, bool Square = (Sub::rows == Sub::columns)>
struct is_symmetric_helper {
  static constexpr bool value = false;
};

template <typename Sub>
struct is_symmetric_helper<Sub, true> {
  template <std::size_t... Ks>
  static constexpr bool all_mirrored(std::index_sequence<Ks...>) {
    constexpr std::size_t C = Sub::columns;
    return (std::same_as<typename Sub::template value_at<Ks / C, Ks % C>,
                         typename Sub::template value_at<Ks % C, Ks / C>> &&
            ...);
  }
  static constexpr bool value =
      all_mirrored(std::make_index_sequence<Sub::rows * Sub::columns>{});
};

template <typename Sub, bool Square = (Sub::rows == Sub::columns)>
struct is_upper_triangular_helper {
  static constexpr bool value = false;
};

template <typename Sub>
struct is_upper_triangular_helper<Sub, true> {
  template <std::size_t... Ks>
  static constexpr bool all_lower_zero(std::index_sequence<Ks...>) {
    constexpr std::size_t C = Sub::columns;
    return (
        (Ks / C <= Ks % C ||
         std::same_as<typename Sub::template value_at<Ks / C, Ks % C>, zero>) &&
        ...);
  }
  static constexpr bool value =
      all_lower_zero(std::make_index_sequence<Sub::rows * Sub::columns>{});
};

template <typename Sub, bool Square = (Sub::rows == Sub::columns)>
struct is_unit_upper_triangular_helper {
  static constexpr bool value = false;
};

template <typename Sub>
struct is_unit_upper_triangular_helper<Sub, true> {
  template <std::size_t... Ks>
  static constexpr bool all_unit_diagonal(std::index_sequence<Ks...>) {
    constexpr std::size_t C = Sub::columns;
    return (
        (Ks / C != Ks % C ||
         std::same_as<typename Sub::template value_at<Ks / C, Ks % C>, one>) &&
        ...);
  }
  static constexpr bool value =
      is_upper_triangular_helper<Sub>::value &&
      all_unit_diagonal(std::make_index_sequence<Sub::rows * Sub::columns>{});
};
}  // namespace detail

template <typename T>
concept is_row = detail::is_row_v<T>;

template <typename T>
concept is_canonical_row = detail::is_canonical_row_v<T>;

template <typename T>
concept is_canonical_substructure = detail::is_canonical_substructure_v<T>;

template <typename T>
concept is_canonical = is_canonical_substructure<T>;

template <typename Sub>
concept is_symmetric =
    is_canonical<Sub> && detail::is_symmetric_helper<Sub>::value;

template <typename Sub>
concept is_upper_triangular =
    is_canonical<Sub> && detail::is_upper_triangular_helper<Sub>::value;

template <typename Sub>
concept is_unit_upper_triangular =
    is_canonical<Sub> && detail::is_unit_upper_triangular_helper<Sub>::value;
}  // namespace checks

namespace getters {
template <typename T>
struct num_values {};

template <typename... value_ts>
struct num_values<row<value_ts...>> {
  static constexpr std::size_t value = sizeof...(value_ts);
};

template <std::size_t N, typename S>
struct substructure_index;

template <std::size_t N, typename... Rs>
struct substructure_index<N, substructure<Rs...>> {
  using type = std::tuple_element_t<N, std::tuple<Rs...>>;
};

template <std::size_t N, typename R>
struct row_index;

template <std::size_t N, typename... Vs>
struct row_index<N, row<Vs...>> {
  using type = std::tuple_element_t<N, std::tuple<Vs...>>;
};

template <typename T, typename...>
struct first {
  using type = T;
};
}  // namespace getters

}  // namespace detail

template <typename... value_ts>
struct row {
  static_assert((detail::checks::is_value<value_ts> && ...));
};

namespace detail {

template <typename V>
constexpr auto negate_value() {
  if constexpr (checks::is_variable<V>) {
    static_assert(checks::is_index_variable<V>);
    return V{};
  } else {
    return integral_value<-V::value>{};
  }
}

template <typename Symbolic, bool FromSameMatrix = false,
          bool CollapseToVariable = true>
constexpr auto resolve_symbolic() {
  static_assert(!(Symbolic::operation == SymbolicOp::DIV &&
                  std::same_as<typename Symbolic::operandB, zero>));

  if constexpr (Symbolic::operation == SymbolicOp::ADD) {
    if constexpr (checks::is_variable<typename Symbolic::operandA> ||
                  checks::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else if constexpr (checks::is_integral_value<
                             typename Symbolic::operandA> &&
                         checks::is_integral_value<
                             typename Symbolic::operandB>) {
      return integral_value<Symbolic::operandA::value +
                            Symbolic::operandB::value>{};
    } else {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::SUB) {
    if constexpr (FromSameMatrix && std::same_as<typename Symbolic::operandA,
                                                 typename Symbolic::operandB>) {
      return zero{};
    } else if constexpr (checks::is_variable<typename Symbolic::operandA> ||
                         checks::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else if constexpr (checks::is_integral_value<
                             typename Symbolic::operandA> &&
                         checks::is_integral_value<
                             typename Symbolic::operandB>) {
      return integral_value<Symbolic::operandA::value -
                            Symbolic::operandB::value>{};
    } else {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::MUL) {
    if constexpr (std::same_as<typename Symbolic::operandA, zero> ||
                  std::same_as<typename Symbolic::operandB, zero>) {
      return zero{};
    } else if constexpr (checks::is_variable<typename Symbolic::operandA> ||
                         checks::is_variable<typename Symbolic::operandB>) {
      if constexpr (CollapseToVariable) {
        return variable{};
      } else {
        return Symbolic{};
      }
    } else {
      return integral_value<Symbolic::operandA::value *
                            Symbolic::operandB::value>{};
    }
  } else if constexpr (Symbolic::operation == SymbolicOp::DIV) {
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

template <typename V, bool = false, bool = true>
struct resolve_value {};

template <typename V, bool FromSameMatrix, bool CollapseToVariable>
  requires(checks::is_symbolic<V>)
struct resolve_value<V, FromSameMatrix, CollapseToVariable> {
  using type =
      decltype(resolve_symbolic<V, FromSameMatrix, CollapseToVariable>());
};
template <typename V, bool FromSameMatrix, bool CollapseToVariable>
  requires(!checks::is_symbolic<V>)
struct resolve_value<V, FromSameMatrix, CollapseToVariable> {
  using type = V;
};

template <typename A, typename B>
  requires(checks::is_value<A> && checks::is_value<B>)
using value_addition = symbolic<SymbolicOp::ADD, A, B>;

template <typename A, typename B>
  requires(checks::is_value<A> && checks::is_value<B>)
using value_subtraction = symbolic<SymbolicOp::SUB, A, B>;

template <typename A, typename B>
  requires(checks::is_value<A> && checks::is_value<B>)
using value_multiplication = symbolic<SymbolicOp::MUL, A, B>;

template <typename A, typename B>
  requires(checks::is_value<A> && checks::is_value<B>)
using value_division = symbolic<SymbolicOp::DIV, A, B>;

template <typename V>
  requires(checks::is_value<V>)
using value_negation = decltype(negate_value<V>());

template <typename T>
struct get_determinant_helper {};

template <typename T>
  requires(T::rows == 1 && T::columns == 1)
struct get_determinant_helper<T> {
  using type = T::template value_at<0, 0>;
};

template <typename T>
  requires(T::rows == 2 && T::columns == 2)
struct get_determinant_helper<T> {
 private:
  using t1 = value_multiplication<typename T::template value_at<0, 0>,
                                  typename T::template value_at<1, 1>>;
  using t2 = value_multiplication<typename T::template value_at<1, 0>,
                                  typename T::template value_at<0, 1>>;

  using r1 = resolve_value<t1, true, false>::type;
  using r2 = resolve_value<t2, true, false>::type;

 public:
  using type = resolve_value<value_subtraction<r1, r2>, true, true>::type;
};

namespace checks {
template <typename T>
concept is_invertible =
    is_canonical<T> && T::rows == T::columns &&
    (T::rows == 1 || T::rows == 2) &&
    !std::same_as<typename get_determinant_helper<T>::type, zero>;
}  // namespace checks

template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(checks::is_canonical<substructure_A> &&
           checks::is_canonical<substructure_B>)
struct get_addition_element {
  using A_value_type = substructure_A::template value_at<I, J>;
  using B_value_type = substructure_B::template value_at<I, J>;
  using type = value_addition<A_value_type, B_value_type>;
};

template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(checks::is_canonical<substructure_A> &&
           checks::is_canonical<substructure_B>)
struct get_subtraction_element {
  using A_value_type = substructure_A::template value_at<I, J>;
  using B_value_type = substructure_B::template value_at<I, J>;
  using type = value_subtraction<A_value_type, B_value_type>;
};

template <typename T, typename S>
struct dot_product {};

template <typename T, typename... Ts, typename S, typename... Ss>
struct dot_product<row<T, Ts...>, row<S, Ss...>> {
  static_assert(sizeof...(Ts) == sizeof...(Ss));

 private:
  using rem_val = typename dot_product<row<Ts...>, row<Ss...>>::type;
  using cur_val = decltype(resolve_symbolic<value_multiplication<T, S>>());

 public:
  using type = decltype(resolve_symbolic<value_addition<cur_val, rem_val>>());
};

template <>
struct dot_product<row<>, row<>> {
  using type = zero;
};

template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(checks::is_canonical<substructure_A> &&
           checks::is_canonical<substructure_B>)
struct get_multiplication_element {
 private:
  using A_row_type = typename substructure_A::template row_at<I>;
  using B_column_type = typename substructure_B::template column_at<J>;

  static_assert(getters::num_values<A_row_type>::value ==
                getters::num_values<B_column_type>::value);

 public:
  using type = typename dot_product<A_row_type, B_column_type>::type;
};

template <typename substructure_A, std::size_t I, std::size_t J>
  requires(checks::is_canonical<substructure_A>)
struct get_negation_element {
  using type = value_negation<typename substructure_A::template value_at<I, J>>;
};

template <typename substructure_A, std::size_t I, std::size_t J>
  requires(checks::is_canonical<substructure_A>)
struct get_identity_element {
  using type = typename substructure_A::template value_at<I, J>;
};

template <template <std::size_t, std::size_t> typename C, std::size_t I,
          typename T>
struct generate_substructure_values {};

template <template <std::size_t, std::size_t> typename C, std::size_t I,
          std::size_t... Js>
struct generate_substructure_values<C, I,
                                    std::integer_sequence<std::size_t, Js...>> {
  using type = row<typename C<I, Js>::type...>;
};

template <template <std::size_t, std::size_t> typename C, std::size_t N,
          typename T>
struct generate_substructure_rows {};

template <template <std::size_t, std::size_t> typename C, std::size_t N,
          std::size_t... Is>
struct generate_substructure_rows<C, N,
                                  std::integer_sequence<std::size_t, Is...>> {
  using type = substructure<typename generate_substructure_values<
      C, Is, std::make_index_sequence<N>>::type...>;
};

namespace builders {

template <std::size_t N>
constexpr std::size_t symmetric_index(std::size_t i, std::size_t j) {
  if (i > j) {
    std::size_t t = i;
    i = j;
    j = t;
  }
  std::size_t idx = 0;
  for (std::size_t r = 0; r < i; ++r)
    idx += (N - r);
  return idx + (j - i);
}

template <std::size_t I, std::size_t J>
struct dense_cell {
  using type = variable;  // every cell an independent variable
};

template <std::size_t I, std::size_t J>
struct diagonal_cell {
  using type = std::conditional_t<(I == J), variable, zero>;
};

template <std::size_t I, std::size_t J>
struct identity_cell {
  using type = std::conditional_t<(I == J), one, zero>;
};

template <std::size_t N>
struct symmetric_cell {
  template <std::size_t I, std::size_t J>
  struct cell {
    using type = index_variable<symmetric_index<N>(I, J)>;
  };
};

template <std::size_t K>
struct left_columns_cell {
  template <std::size_t I, std::size_t J>
  struct cell {
    using type = std::conditional_t<(J < K), variable, zero>;
  };
};

}  // namespace builders

}  // namespace detail

template <template <std::size_t, std::size_t> typename C, std::size_t Rows,
          std::size_t Cols>
using make_substructure = typename detail::generate_substructure_rows<
    C, Cols, std::make_integer_sequence<std::size_t, Rows>>::type;

// Rows x Cols, every cell an independent variable (structure-blind dense).
template <std::size_t Rows, std::size_t Cols = Rows>
using make_dense_substructure =
    make_substructure<detail::builders::dense_cell, Rows, Cols>;

// NxN symmetric: (i, j) and (j, i) name the same variable.
template <std::size_t N>
using make_symmetric_substructure =
    make_substructure<detail::builders::symmetric_cell<N>::template cell, N, N>;

// NxN diagonal: variables on the diagonal, structural zeros off it.
template <std::size_t N>
using make_diagonal_substructure =
    make_substructure<detail::builders::diagonal_cell, N, N>;

// NxN identity: structural ones on the diagonal, zeros off it.
template <std::size_t N>
using make_identity_substructure =
    make_substructure<detail::builders::identity_cell, N, N>;

// Rows x Cols whose leftmost K columns are variables and the rest zero.
template <std::size_t Rows, std::size_t Cols, std::size_t K>
using make_left_columns_substructure =
    make_substructure<detail::builders::left_columns_cell<K>::template cell,
                      Rows, Cols>;

namespace detail {

template <typename A, typename B>
  requires(checks::is_canonical<A> && checks::is_canonical<B>)
struct additive_substructure {
  static_assert(A::columns == B::columns);
  static_assert(A::rows == B::rows);

  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = get_addition_element<A, B, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename A, typename B>
  requires(checks::is_canonical<A> && checks::is_canonical<B>)
struct subtractive_substructure {
  static_assert(A::columns == B::columns);
  static_assert(A::rows == B::rows);

  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = get_subtraction_element<A, B, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename A, typename B>
  requires(checks::is_canonical<A> && checks::is_canonical<B>)
struct multiplicative_substructure {
  static_assert(A::columns == B::rows);

  static constexpr std::size_t columns = B::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = get_multiplication_element<A, B, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename A>
  requires(checks::is_canonical<A>)
struct negated_substructure {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = get_negation_element<A, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename A>
  requires((A::columns == A::rows))
struct upper_symmetric_substructure {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = std::conditional_t<(I < J), get_identity_element<A, I, J>,
                                      get_identity_element<A, J, I>>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename T, bool Sym>
struct upper_symmetric_if {};

template <typename T>
struct upper_symmetric_if<T, false> {
  using type = T;
};

template <typename T>
struct upper_symmetric_if<T, true> {
  using type = typename upper_symmetric_substructure<T>::type;
};

template <typename A>
struct inverted_substructure {};

template <typename A>
  requires(checks::is_canonical<A> && checks::is_invertible<A> && A::rows == 1)
struct inverted_substructure<A> {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  struct accessor {
    using type = value_division<integral_value<1>,
                                typename get_identity_element<A, 0, 0>::type>;
  };

  using raw_type = typename generate_substructure_rows<
      accessor, columns, std::make_index_sequence<rows>>::type::canonical_type;

  using type = upper_symmetric_if<raw_type, checks::is_symmetric<A>>::type;
};

template <typename A>
  requires(checks::is_canonical<A> && checks::is_invertible<A> && A::rows == 2)
struct inverted_substructure<A> {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  using det = get_determinant_helper<A>::type;

  template <std::size_t I, std::size_t J>
  struct accessor {
    using type = decltype([]() {
      if constexpr (I == 0 && J == 0) {
        return value_division<typename get_identity_element<A, 1, 1>::type,
                              det>{};
      } else if constexpr (I == 1 && J == 1) {
        return value_division<typename get_identity_element<A, 0, 0>::type,
                              det>{};
      } else {
        using neg = typename resolve_value<
            value_subtraction<zero,
                              typename get_identity_element<A, I, J>::type>,
            false, false>::type;
        return value_division<neg, det>{};
      }
    }());
  };

  using raw_type = typename generate_substructure_rows<
      accessor, columns, std::make_index_sequence<rows>>::type::canonical_type;

  using type = upper_symmetric_if<raw_type, checks::is_symmetric<A>>::type;
};

template <typename A>
  requires(checks::is_canonical<A>)
struct transposed_substructure {
  static constexpr std::size_t columns = A::rows;
  static constexpr std::size_t rows = A::columns;

  template <std::size_t I, std::size_t J>
  using accessor = get_identity_element<A, J, I>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

template <typename A, typename B>
  requires(checks::is_canonical<A> && checks::is_canonical<B>)
struct congruence_substructure {
  static_assert(A::columns == B::rows);
  static_assert(B::columns == A::columns);

  using AB_type =
      typename multiplicative_substructure<A, B>::type::canonical_type;
  using ABAT_type = typename multiplicative_substructure<
      AB_type, typename transposed_substructure<A>::type::canonical_type>::
      type::canonical_type;

  using sym_ABAT_type =
      typename upper_symmetric_substructure<ABAT_type>::type::canonical_type;
  static_assert(checks::is_symmetric<sym_ABAT_type>);

  static constexpr std::size_t columns = A::rows;
  static constexpr std::size_t rows = A::rows;

  using type = sym_ABAT_type;
};

template <typename substructure_t, std::size_t J, typename T>
struct column_getter_helper {};

template <typename substructure_t, std::size_t J, std::size_t... Is>
struct column_getter_helper<substructure_t, J,
                            std::integer_sequence<std::size_t, Is...>> {
  using type = row<typename substructure_t::template value_at<Is, J>...>;
};

template <typename T>
struct index_variable_getter {
  using type = std::integral_constant<long, -1>;
};

template <std::size_t I>
struct index_variable_getter<index_variable<I>> {
  using type = std::integral_constant<long, I>;
};

template <std::size_t N>
struct index_table_t {
  std::array<std::size_t, N> table;
  std::size_t num_variables;
};

template <typename Sub, std::size_t... Ks>
  requires(checks::is_canonical<Sub>)
constexpr auto make_index_table(std::index_sequence<Ks...>) {
  constexpr std::size_t C = Sub::columns;
  constexpr std::array<bool, sizeof...(Ks)> is_indexvar = {
      (checks::is_index_variable<
          typename Sub::template value_at<Ks / C, Ks % C>>)...};
  constexpr std::array<long, sizeof...(Ks)> indexvar_indices = {
      (index_variable_getter<
          typename Sub::template value_at<Ks / C, Ks % C>>::type::value)...};

  std::array<std::size_t, sizeof...(Ks)> table{};
  std::size_t acc = 0;
  for (std::size_t k = 0; k < sizeof...(Ks); ++k) {
    if (is_indexvar[k]) {
      bool seen_prev = false;
      std::size_t prev_idx = 0;
      for (std::size_t j = 0; j < k; ++j) {
        if (is_indexvar[j] && indexvar_indices[k] == indexvar_indices[j]) {
          seen_prev = true;
          prev_idx = j;
          break;
        }
      }

      if (!seen_prev) {
        table[k] = acc;
        ++acc;
      } else {
        table[k] = table[prev_idx];
      }
    } else {
      table[k] = std::numeric_limits<std::size_t>::max();
    }
  }
  return index_table_t<sizeof...(Ks)>{table, acc};
}

template <typename Sub>
  requires(checks::is_canonical<Sub>)
inline constexpr auto index_table_info =
    make_index_table<Sub>(std::make_index_sequence<Sub::rows * Sub::columns>{});

template <typename Sub>
  requires(checks::is_canonical<Sub>)
inline constexpr auto index_table = index_table_info<Sub>.table;

template <typename Sub>
  requires(checks::is_canonical<Sub>)
inline constexpr auto get_num_variables = index_table_info<Sub>.num_variables;

template <typename Sub, typename Val, std::size_t... Ks>
constexpr std::size_t first_equal_flat_index(std::index_sequence<Ks...>) {
  constexpr std::size_t C = Sub::columns;
  constexpr std::array<bool, sizeof...(Ks)> eq = {
      std::same_as<typename Sub::template value_at<Ks / C, Ks % C>, Val>...};
  for (std::size_t k = 0; k < sizeof...(Ks); ++k) {
    if (eq[k]) {
      return k;
    }
  }
  return sizeof...(Ks);
}

template <typename Sub, typename Val>
inline constexpr std::size_t first_equal_flat_index_v =
    first_equal_flat_index<Sub, Val>(
        std::make_index_sequence<Sub::rows * Sub::columns>{});

template <typename T, std::size_t I, std::size_t J>
struct canonicalize_substructure_functor {
 private:
  using curr_value = T::template value_at<I, J>;
  using resolved_value = typename resolve_value<curr_value>::type;

 public:
  using type = std::conditional_t<
      std::same_as<curr_value, variable>, index_variable<T::columns * I + J>,
      std::conditional_t<
          checks::is_integral_value<resolved_value>, resolved_value,
          index_variable<first_equal_flat_index_v<T, curr_value>>>>;
};

template <typename T>
struct canonicalize_substructure {
  template <std::size_t I, std::size_t J>
  using accessor = canonicalize_substructure_functor<T, I, J>;

  using type = typename generate_substructure_rows<
      accessor, T::columns, std::make_index_sequence<T::rows>>::type;
};

}  // namespace detail

template <typename... row_ts>
struct substructure {
  using this_t = substructure<row_ts...>;

  static constexpr std::size_t columns = detail::getters::num_values<
      typename detail::getters::first<row_ts...>::type>::value;
  static constexpr std::size_t rows = sizeof...(row_ts);
  static_assert(((detail::getters::num_values<row_ts>::value == columns) &&
                 ...),
                "Rows do not all have the same size!");

  template <std::size_t I, std::size_t J>
  using value_at = typename detail::getters::row_index<
      J, typename detail::getters::substructure_index<
             I, substructure<row_ts...>>::type>::type;

  using canonical_type =
      typename detail::canonicalize_substructure<this_t>::type;

  template <std::size_t I, std::size_t J>
  static constexpr std::size_t element_index =
      detail::index_table<canonical_type>[I * columns + J];
  static constexpr std::size_t num_variables =
      detail::get_num_variables<canonical_type>;

  static_assert((detail::checks::is_row<row_ts> && ...));

  template <std::size_t I>
  using row_at = typename detail::getters::substructure_index<I, this_t>::type;

  template <std::size_t J>
  using column_at = typename detail::column_getter_helper<
      substructure<row_ts...>, J, std::make_index_sequence<rows>>::type;

  template <typename other_substructure_t>
    requires(detail::checks::is_canonical<this_t> &&
             detail::checks::is_canonical<other_substructure_t>)
  using addition_type =
      detail::additive_substructure<this_t,
                                    other_substructure_t>::type::canonical_type;

  template <typename other_substructure_t>
    requires(detail::checks::is_canonical<this_t> &&
             detail::checks::is_canonical<other_substructure_t>)
  using subtraction_type = detail::subtractive_substructure<
      this_t, other_substructure_t>::type::canonical_type;

  template <typename other_substructure_t>
    requires(detail::checks::is_canonical<this_t> &&
             detail::checks::is_canonical<other_substructure_t>)
  using multiplication_type = detail::multiplicative_substructure<
      this_t, other_substructure_t>::type::canonical_type;

  template <typename other_substructure_t>
    requires(detail::checks::is_canonical<this_t> &&
             detail::checks::is_canonical<other_substructure_t> &&
             detail::checks::is_symmetric<other_substructure_t>)
  using congruence_type = detail::congruence_substructure<
      this_t, other_substructure_t>::type::canonical_type;
};

template <typename substructure_t, typename scalar_t>
struct matrix {
  using substructure_type = substructure_t;
  using canonical_substructure_type = typename substructure_t::canonical_type;

  static_assert(substructure_type::columns ==
                canonical_substructure_type::columns);
  static_assert(substructure_type::rows == canonical_substructure_type::rows);
  static_assert(std::same_as<canonical_substructure_type,
                             typename detail::canonicalize_substructure<
                                 canonical_substructure_type>::type>);

  static_assert(detail::checks::is_canonical<canonical_substructure_type>);

  using scalar_type = scalar_t;

  static constexpr std::size_t rows = substructure_t::rows;
  static constexpr std::size_t columns = substructure_t::columns;

  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_t at() const {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    static_assert(detail::checks::is_canonical_value<value_t>);

    if constexpr (detail::checks::is_index_variable<value_t>) {
      return m_elems[canonical_substructure_type::template element_index<I, J>];
    } else {
      return static_cast<scalar_t>(
          canonical_substructure_type::template value_at<I, J>::value);
    }
  }

  template <std::size_t I, std::size_t J>
    requires(detail::checks::is_index_variable<
             typename canonical_substructure_type::template value_at<I, J>>)
  DETRAY_HOST_DEVICE scalar_t &at() {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    static_assert(detail::checks::is_canonical_value<value_t>);
    return m_elems[canonical_substructure_type::template element_index<I, J>];
  }

  /// @name Compile-time element access in detray's spelling
  ///
  /// detray's generated code reaches into matrices through
  /// @c getter::element<I, J>, which dispatches to a member @c element<I, J>
  /// when the type has one. Providing it lets a @c ksm::matrix be handed to
  /// that code unchanged. These are only another spelling of @c at.
  /// @{
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_t element() const {
    return this->template at<I, J>();
  }

  template <std::size_t I, std::size_t J>
    requires(detail::checks::is_index_variable<
             typename canonical_substructure_type::template value_at<I, J>>)
  DETRAY_HOST_DEVICE scalar_t &element() {
    return this->template at<I, J>();
  }
  /// @}

  /// @returns the identity matrix in this substructure
  ///
  /// @c matrix::identity cannot serve here: it writes the diagonal through
  /// runtime indices, and a structured matrix is only addressable at compile
  /// time. Well formed only when the substructure can represent the identity
  /// at all, i.e. when every structural cell already holds the value the
  /// identity requires there. That is checked at compile time.
  DETRAY_HOST_DEVICE static matrix identity() {
    matrix rv{};
    rv.set_identity_cells(std::make_index_sequence<rows * columns>{});
    return rv;
  }

  /// @returns this matrix as a dense detray matrix of the algebra
  /// @tparam algebra_t
  ///
  /// Structural cells materialise their compile-time value, so the result is a
  /// full @c rows x @c columns matrix with nothing left implicit. This is the
  /// escape hatch that lets a structured matrix be handed to code that has not
  /// been converted, so that call sites can be migrated one at a time.
  ///
  /// @note Convert at the edges of a computation, never between the operations
  /// inside one. A round trip through a dense matrix discards exactly the
  /// structure that would have made the next operation cheap.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dmatrix<algebra_t, rows, columns> to_dense() const {
    dmatrix<algebra_t, rows, columns> rv{};
    to_dense_cells<algebra_t>(rv, std::make_index_sequence<rows * columns>{});
    return rv;
  }

  /// @returns a matrix of this substructure carrying the values of the dense
  /// matrix @p m
  ///
  /// Only the free cells are read: a structural cell has no storage, so a
  /// dense matrix that does not actually have the claimed structure cannot be
  /// represented, and whatever those cells held would be silently discarded.
  /// Debug builds therefore check the claim at the one place it can be
  /// checked -- this boundary. Every constant cell must hold its compile-time
  /// value, and cells that share a variable must agree with one another. This
  /// is the same discipline the sympy-generated
  /// @c transport_covariance_to_bound_impl already applies to the full
  /// Jacobian, where it opens with eleven such assertions.
  ///
  /// @note See @c to_dense on keeping conversions at the edges of a
  /// computation rather than between the operations inside one.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE static matrix from_dense(
      const dmatrix<algebra_t, rows, columns> &m) {
    matrix rv{};
    rv.from_dense_cells(m, std::make_index_sequence<rows * columns>{});
    return rv;
  }

  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator*(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type = matrix<
        typename canonical_substructure_type::template multiplication_type<
            typename other_substructure_t::canonical_type>,
        scalar_t>;

    result_type rv{};
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      return this->template product_dot<I, J>(
          o, std::make_index_sequence<canonical_substructure_type::columns>{});
    });
    return rv;
  }

  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator+(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template addition_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    result_type rv{};
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      using A_val =
          typename canonical_substructure_type::template value_at<I, J>;
      using B_val =
          typename other_substructure_t::canonical_type::template value_at<I,
                                                                           J>;
      if constexpr (std::same_as<A_val, zero>) {
        return o.template at<I, J>();
      } else if constexpr (std::same_as<B_val, zero>) {
        return this->template at<I, J>();
      } else {
        return this->template at<I, J>() + o.template at<I, J>();
      }
    });
    return rv;
  }

  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator-(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template subtraction_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    result_type rv{};
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      using A_val =
          typename canonical_substructure_type::template value_at<I, J>;
      using B_val =
          typename other_substructure_t::canonical_type::template value_at<I,
                                                                           J>;
      if constexpr (std::same_as<A_val, zero>) {
        return -o.template at<I, J>();
      } else if constexpr (std::same_as<B_val, zero>) {
        return this->template at<I, J>();
      } else {
        return this->template at<I, J>() - o.template at<I, J>();
      }
    });
    return rv;
  }

  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto congruence(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template congruence_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    const auto ab = *this * o;
    const auto at = this->transpose();

    result_type rv{};
    fill_with(rv, [this, &o, &ab, &at]<std::size_t I, std::size_t J>() {
      return ab.template product_dot<I, J>(
          at, std::make_index_sequence<canonical_substructure_type::columns>{});
    });
    return rv;
  }

  DETRAY_HOST_DEVICE auto inverse() const {
    using result_type = matrix<typename detail::inverted_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};

    static_assert(rows == 1 || rows == 2);
    static_assert(detail::checks::is_invertible<canonical_substructure_type>);

    using determinant_type =
        detail::get_determinant_helper<canonical_substructure_type>::type;

    static_assert(!std::same_as<determinant_type, zero>);

    scalar_t det;

    if constexpr (std::same_as<determinant_type, one>) {
      // This should never be used.
      det = scalar_t{0};
    } else {
      if constexpr (rows == 1) {
        det = this->template at<0, 0>();
      } else {
        det = this->template at<0, 0>() * this->template at<1, 1>() -
              this->template at<1, 0>() * this->template at<0, 1>();
      }
    }

    if constexpr (rows == 1) {
      fill_with(rv, [this, det]<std::size_t I, std::size_t J>() {
        if constexpr (std::same_as<determinant_type, one>) {
          return this->template at<I, J>();
        } else {
          return scalar_t{1} / det;
        }
      });
    } else {
      fill_with(rv, [this, det]<std::size_t I, std::size_t J>() {
        if constexpr (std::same_as<determinant_type, one>) {
          if constexpr (I == 0 && J == 0) {
            return this->template at<1, 1>();
          } else if constexpr (I == 1 && J == 1) {
            return this->template at<0, 0>();
          } else {
            return -this->template at<I, J>();
          }
        } else {
          if constexpr (I == 0 && J == 0) {
            return this->template at<1, 1>() / det;
          } else if constexpr (I == 1 && J == 1) {
            return this->template at<0, 0>() / det;
          } else {
            return -this->template at<I, J>() / det;
          }
        }
      });
    }

    return rv;
  }

  DETRAY_HOST_DEVICE auto transpose() const {
    using result_type = matrix<typename detail::transposed_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};
    fill_with(rv, [this]<std::size_t I, std::size_t J>() {
      return this->template at<J, I>();
    });
    return rv;
  }

  DETRAY_HOST_DEVICE auto operator-() const {
    using result_type = matrix<typename detail::negated_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};
    fill_with(rv, [this]<std::size_t I, std::size_t J>() {
      return -this->template at<I, J>();
    });
    return rv;
  }

 private:
  /// Give cell (I, J) the value the identity matrix has there, or check that
  /// the substructure already fixes it to that value.
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE void set_identity_cell() {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    if constexpr (detail::checks::is_index_variable<value_t>) {
      constexpr std::size_t first_occurrence =
          detail::first_equal_flat_index_v<canonical_substructure_type,
                                           value_t>;

      if constexpr (first_occurrence == I * columns + J) {
        this->template at<I, J>() = (I == J) ? scalar_t{1} : scalar_t{0};
      } else {
        // Shares storage with an earlier cell. The identity is representable
        // only if the two want the same value there, which a substructure
        // pairing a diagonal cell with an off-diagonal one would not.
        assert((this->template at<I, J>() ==
                ((I == J) ? scalar_t{1} : scalar_t{0})));
      }
    } else {
      static_assert(value_t::value == ((I == J) ? 1 : 0),
                    "this substructure cannot represent the identity matrix: a "
                    "structural cell fixes a value the identity does not have");
    }
  }

  template <std::size_t... Ks>
  DETRAY_HOST_DEVICE void set_identity_cells(std::index_sequence<Ks...>) {
    (set_identity_cell<Ks / columns, Ks % columns>(), ...);
  }

  /// Write every cell, structural ones included, into the dense matrix @p rv.
  /// The scalar cast is explicit because the algebra's scalar need not be the
  /// scalar this matrix stores.
  template <typename algebra_t, typename dense_t, std::size_t... Ks>
  DETRAY_HOST_DEVICE void to_dense_cells(dense_t &rv,
                                         std::index_sequence<Ks...>) const {
    ((getter::element<Ks / columns, Ks % columns>(rv) =
          static_cast<dscalar<algebra_t>>(
              this->template at<Ks / columns, Ks % columns>())),
     ...);
  }

  /// Take cell (I, J) of the dense matrix @p m: store it if the substructure
  /// declares it a free value, otherwise check it against what the
  /// substructure promises, since there is nowhere to put a differing value.
  template <std::size_t I, std::size_t J, typename dense_t>
  DETRAY_HOST_DEVICE void from_dense_cell([[maybe_unused]] const dense_t &m) {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    if constexpr (detail::checks::is_index_variable<value_t>) {
      constexpr std::size_t first_occurrence =
          detail::first_equal_flat_index_v<canonical_substructure_type,
                                           value_t>;

      if constexpr (first_occurrence == I * columns + J) {
        this->template at<I, J>() =
            static_cast<scalar_t>(getter::element<I, J>(m));
      } else {
        // Shares storage with an earlier cell. The substructure promises the
        // two are the same value -- a dense input where they differ is not the
        // matrix it claims to be, and one of the two would be dropped.
        // (Doubled parentheses: the commas in the template argument lists
        // would otherwise be read as macro argument separators.)
        assert((static_cast<scalar_t>(getter::element<I, J>(m)) ==
                this->template at<I, J>()));
      }
    } else {
      assert((static_cast<scalar_t>(getter::element<I, J>(m)) ==
              static_cast<scalar_t>(value_t::value)));
    }
  }

  /// Cells are visited in flat order, so the first occurrence of a shared
  /// variable is always stored before the cells that have to agree with it.
  template <typename dense_t, std::size_t... Ks>
  DETRAY_HOST_DEVICE void from_dense_cells(const dense_t &m,
                                           std::index_sequence<Ks...>) {
    (from_dense_cell<Ks / columns, Ks % columns>(m), ...);
  }

  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t K>
  DETRAY_HOST_DEVICE static constexpr long dot_term_constant() {
    using A_val = typename canonical_substructure_type::template value_at<I, K>;
    using B_val = typename OtherSub::canonical_type::template value_at<K, J>;
    if constexpr (detail::checks::is_integral_value<A_val> &&
                  detail::checks::is_integral_value<B_val>) {
      return static_cast<long>(A_val::value) * B_val::value;
    } else {
      return 0;
    }
  }

  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t... Ks>
  DETRAY_HOST_DEVICE static constexpr long dot_constant_offset(
      std::index_sequence<Ks...>) {
    return (0 + ... + dot_term_constant<I, J, OtherSub, Ks>());
  }

  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t K,
            std::size_t... Ks>
  DETRAY_HOST_DEVICE auto product_dot_vars(
      const matrix<OtherSub, scalar_t> &o,
      std::index_sequence<K, Ks...>) const {
    using A_val = typename canonical_substructure_type::template value_at<I, K>;
    using B_val = typename OtherSub::canonical_type::template value_at<K, J>;

    constexpr bool is_zero_term =
        std::same_as<A_val, zero> || std::same_as<B_val, zero>;
    constexpr bool both_const = detail::checks::is_integral_value<A_val> &&
                                detail::checks::is_integral_value<B_val>;
    constexpr bool is_var_term = !is_zero_term && !both_const;

    if constexpr (sizeof...(Ks) > 0) {
      const auto rem = product_dot_vars<I, J>(o, std::index_sequence<Ks...>{});
      using rem_type = std::decay_t<decltype(rem)>;

      if constexpr (!is_var_term) {
        return rem;
      } else if constexpr (std::same_as<A_val, one>) {
        if constexpr (std::same_as<rem_type, zero>) {
          return o.template at<K, J>();
        } else {
          return o.template at<K, J>() + rem;
        }
      } else if constexpr (std::same_as<B_val, one>) {
        if constexpr (std::same_as<rem_type, zero>) {
          return this->template at<I, K>();
        } else {
          return this->template at<I, K>() + rem;
        }
      } else {
        if constexpr (std::same_as<rem_type, zero>) {
          return this->template at<I, K>() * o.template at<K, J>();
        } else {
          return this->template at<I, K>() * o.template at<K, J>() + rem;
        }
      }
    } else {
      if constexpr (!is_var_term) {
        return zero{};
      } else if constexpr (std::same_as<A_val, one>) {
        return o.template at<K, J>();
      } else if constexpr (std::same_as<B_val, one>) {
        return this->template at<I, K>();
      } else {
        return this->template at<I, K>() * o.template at<K, J>();
      }
    }
  }

  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t... Ks>
  DETRAY_HOST_DEVICE auto product_dot(const matrix<OtherSub, scalar_t> &o,
                                      std::index_sequence<Ks...> seq) const {
    constexpr long offset = dot_constant_offset<I, J, OtherSub>(seq);
    const auto var_sum = product_dot_vars<I, J>(o, seq);
    using vs_type = std::decay_t<decltype(var_sum)>;

    if constexpr (std::same_as<vs_type, zero>) {
      return scalar_t(offset);
    } else if constexpr (offset == 0) {
      return var_sum;
    } else {
      return var_sum + scalar_t(offset);
    }
  }

  template <std::size_t I, std::size_t J, typename Result, typename F>
  DETRAY_HOST_DEVICE void fill_cell(Result &rv, F &f) const {
    if constexpr (requires {
                    rv.template at<I, J>() = f.template operator()<I, J>();
                  }) {
      if constexpr (detail::first_equal_flat_index_v<
                        typename Result::canonical_substructure_type,
                        typename Result::canonical_substructure_type::
                            template value_at<I, J>> ==
                    I * Result::columns + J) {
        rv.template at<I, J>() = f.template operator()<I, J>();
      }
    }
  }

  template <typename Result, typename F, std::size_t... Ks>
  DETRAY_HOST_DEVICE void fill_all(Result &rv, F f,
                                   std::index_sequence<Ks...>) const {
    constexpr std::size_t C = Result::canonical_substructure_type::columns;
    (fill_cell<Ks / C, Ks % C>(rv, f), ...);
  }

  template <typename Result, typename F>
  DETRAY_HOST_DEVICE void fill_with(Result &rv, F f) const {
    using result_sub = typename Result::canonical_substructure_type;
    fill_all(
        rv, f,
        std::make_index_sequence<result_sub::rows * result_sub::columns>{});
  }

  std::array<scalar_t, canonical_substructure_type::num_variables> m_elems{};

  template <class, class>
  friend class matrix;
};

namespace detail {

/// Read cell (I, J) of @p src, which may be a dense matrix or a vector. A
/// vector stands for a single column, so only its row index is meaningful.
template <std::size_t I, std::size_t J, typename src_t>
DETRAY_HOST_DEVICE decltype(auto) element_of(const src_t &src) {
  if constexpr (::detray::concepts::matrix<src_t>) {
    return ::detray::getter::element<I, J>(src);
  } else {
    static_assert(J == 0u, "a vector stands for a single column");
    return ::detray::getter::element<I>(src);
  }
}

/// Copy one element of @p src into cell (DstRow, DstCol) of @p dst.
///
/// A copied region will in general overlap cells the substructure declares
/// structural -- the direction/angle blocks of the bound-to-free Jacobian, for
/// instance, have a zero in the corner. Such a cell has no storage, so the
/// value cannot be written; instead the source has to agree with what the
/// substructure already fixes there, which is checked in debug builds. This is
/// the same discipline @c matrix::from_dense applies at its boundary.
template <std::size_t DstRow, std::size_t DstCol, std::size_t SrcRow,
          std::size_t SrcCol, typename dst_t, typename src_t>
DETRAY_HOST_DEVICE void copy_one_element(dst_t &dst, const src_t &src) {
  using value_t =
      typename dst_t::canonical_substructure_type::template value_at<DstRow,
                                                                     DstCol>;
  using scalar_t = typename dst_t::scalar_type;

  [[maybe_unused]] const auto value =
      static_cast<scalar_t>(element_of<SrcRow, SrcCol>(src));

  if constexpr (checks::is_index_variable<value_t>) {
    dst.template at<DstRow, DstCol>() = value;
  } else {
    assert((value == static_cast<scalar_t>(value_t::value)));
  }
}

}  // namespace detail

/// Copy every element of @p src into @p dst, placing the top-left corner of
/// the copy at (RowOffset, ColOffset).
///
/// This is deliberately not a block assignment. A block would have to be a
/// view over storage that a substructure need not lay out contiguously, and
/// need not keep at all. Copying element by element sidesteps that: free cells
/// are written, and cells the substructure fixes are checked against the
/// source rather than silently overwritten or discarded.
///
/// @p src may be a dense matrix or a vector; a vector is taken as one column.
template <std::size_t RowOffset, std::size_t ColOffset, typename substructure_t,
          typename scalar_t, typename src_t>
DETRAY_HOST_DEVICE void copy_elements_into(
    matrix<substructure_t, scalar_t> &dst, const src_t &src) {
  constexpr auto src_rows =
      static_cast<std::size_t>(::detray::traits::rows<src_t>);
  constexpr auto src_columns =
      static_cast<std::size_t>(::detray::traits::columns<src_t>);

  static_assert(RowOffset + src_rows <= substructure_t::rows,
                "the copied elements do not fit the destination");
  static_assert(ColOffset + src_columns <= substructure_t::columns,
                "the copied elements do not fit the destination");

  [&dst, &src]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (detail::copy_one_element<RowOffset + Ks / src_columns,
                              ColOffset + Ks % src_columns, Ks / src_columns,
                              Ks % src_columns>(dst, src),
     ...);
  }(std::make_index_sequence<src_rows * src_columns>{});
}

}  // namespace detray::ksm

namespace detray::traits {

/// Teach detray's matrix traits about @c ksm::matrix, so that it satisfies
/// @c concepts::matrix and @c concepts::square_matrix and can be passed to
/// generated code whose constraints are written in terms of
/// @c traits::rows / @c columns / @c max_rank. Together with the member
/// @c element<I, J>, this is everything such code needs.
/// @{
template <typename substructure_t, typename scalar_t>
struct index<::detray::ksm::matrix<substructure_t, scalar_t>> {
  using type = std::size_t;
};

template <typename substructure_t, typename scalar_t>
struct dimensions<::detray::ksm::matrix<substructure_t, scalar_t>> {
  using size_type = std::size_t;

  static constexpr size_type _dim{2};
  static constexpr size_type _rows{substructure_t::rows};
  static constexpr size_type _columns{substructure_t::columns};
};
/// @}

}  // namespace detray::traits
