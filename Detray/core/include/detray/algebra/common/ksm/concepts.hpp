// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/detail/symbolic_definitions.hpp"
#include "detray/algebra/common/ksm/fwd.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <concepts>
#include <cstddef>
#include <utility>

// This file defines many of the concepts KSM supplies, and some helper
// traits.
namespace detray::ksm {
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

namespace concepts {
/// @brief True iff the value is a constant integral value.
template <typename T>
concept is_integral_value = detail::is_integral_value_helper<T>::value;

/// @brief True iff the value is a symbolic operation.
template <typename T>
concept is_symbolic = detail::is_symbolic_helper<T>::value;

/// @brief True iff the value is an indexed (i.e. non-anonymous) variable.
template <typename T>
concept is_index_variable = detail::is_index_variable_helper<T>::value;

/// @brief True iff the value is an anonymous or non-anonymous variable.
template <typename T>
concept is_variable = is_index_variable<T> || std::same_as<T, variable>;

/// @brief True iff the type is any valid KSM value, i.e. a symbolic
/// operation, an anonymous variable, a non-anonymous variable, or an integral
/// constant.
template <typename T>
concept is_value = is_symbolic<T> || is_integral_value<T> || is_variable<T>;

/// @brief True iff the type is a canonical KSM value, i.e. an integral
/// constant or a non-anonymous variable.
template <typename T>
concept is_canonical_value = is_integral_value<T> || is_index_variable<T>;
}  // namespace concepts

namespace detail {
template <typename T>
inline constexpr bool is_row_v = false;
template <typename... Vs>
inline constexpr bool is_row_v<row<Vs...>> = (concepts::is_value<Vs> && ...);

template <typename T>
inline constexpr bool is_canonical_row_v = false;
template <typename... Vs>
inline constexpr bool is_canonical_row_v<row<Vs...>> =
    (concepts::is_canonical_value<Vs> && ...);

template <typename T>
inline constexpr bool is_canonical_substructure_v = false;
template <typename... Rs>
inline constexpr bool is_canonical_substructure_v<substructure<Rs...>> =
    (is_canonical_row_v<Rs> && ...);

template <typename substructure_t,
          bool Square = (substructure_t::rows == substructure_t::columns)>
struct is_symmetric_helper {
  static constexpr bool value = false;
};

template <typename substructure_t>
struct is_symmetric_helper<substructure_t, true> {
  template <std::size_t... Ks>
  static constexpr bool all_mirrored(std::index_sequence<Ks...>) {
    constexpr std::size_t C = substructure_t::columns;
    return (std::same_as<
                typename substructure_t::template value_at<Ks / C, Ks % C>,
                typename substructure_t::template value_at<Ks % C, Ks / C>> &&
            ...);
  }
  static constexpr bool value =
      all_mirrored(std::make_index_sequence<substructure_t::rows *
                                            substructure_t::columns>{});
};

template <typename substructure_t,
          bool Square = (substructure_t::rows == substructure_t::columns)>
struct is_diagonal_helper {
  static constexpr bool value = false;
};

template <typename substructure_t>
struct is_diagonal_helper<substructure_t, true> {
  template <std::size_t... Ks>
  static constexpr bool all_off_diagonal_zero(std::index_sequence<Ks...>) {
    constexpr std::size_t C = substructure_t::columns;
    return ((Ks / C == Ks % C ||
             std::same_as<
                 typename substructure_t::template value_at<Ks / C, Ks % C>,
                 zero>) &&
            ...);
  }
  static constexpr bool value = all_off_diagonal_zero(
      std::make_index_sequence<substructure_t::rows *
                               substructure_t::columns>{});
};

template <typename substructure_t,
          bool Square = (substructure_t::rows == substructure_t::columns)>
struct is_upper_triangular_helper {
  static constexpr bool value = false;
};

template <typename substructure_t>
struct is_upper_triangular_helper<substructure_t, true> {
  template <std::size_t... Ks>
  static constexpr bool all_lower_zero(std::index_sequence<Ks...>) {
    constexpr std::size_t C = substructure_t::columns;
    return ((Ks / C <= Ks % C ||
             std::same_as<
                 typename substructure_t::template value_at<Ks / C, Ks % C>,
                 zero>) &&
            ...);
  }
  static constexpr bool value =
      all_lower_zero(std::make_index_sequence<substructure_t::rows *
                                              substructure_t::columns>{});
};

template <typename substructure_t,
          bool Square = (substructure_t::rows == substructure_t::columns)>
struct is_unit_upper_triangular_helper {
  static constexpr bool value = false;
};

/// @brief True iff cell (I, J) already holds what the identity needs there.
///
/// A stored cell can simply be assigned, so only a constant one can fail.
///
/// @{
template <typename substructure_t, std::size_t I, std::size_t J,
          bool Stored = concepts::is_index_variable<
              typename substructure_t::template value_at<I, J>>>
struct identity_constant_ok {
  static constexpr bool value =
      substructure_t::template value_at<I, J>::value == ((I == J) ? 1 : 0);
};

template <typename substructure_t, std::size_t I, std::size_t J>
struct identity_constant_ok<substructure_t, I, J, true> {
  static constexpr bool value = true;
};
/// @}

/// @brief True iff flat cells A and B do not disagree about the identity.
///
/// Two cells naming the same variable share one scalar, so they cannot hold
/// different values. That rules out a substructure pairing a diagonal cell
/// with an off-diagonal one, since the identity wants one in the first and
/// zero in the second.
template <typename substructure_t, std::size_t A, std::size_t B>
struct identity_sharing_ok {
 private:
  static constexpr std::size_t C = substructure_t::columns;
  using a_value = typename substructure_t::template value_at<A / C, A % C>;
  using b_value = typename substructure_t::template value_at<B / C, B % C>;

 public:
  static constexpr bool value = !(concepts::is_index_variable<a_value> &&
                                  std::same_as<a_value, b_value>) ||
                                ((A / C == A % C) == (B / C == B % C));
};

template <typename substructure_t>
struct can_represent_identity_helper {
  static constexpr std::size_t cells =
      substructure_t::rows * substructure_t::columns;

  template <std::size_t... Ks>
  static constexpr bool all_constants_ok(std::index_sequence<Ks...>) {
    constexpr std::size_t C = substructure_t::columns;
    return (identity_constant_ok<substructure_t, Ks / C, Ks % C>::value && ...);
  }

  template <std::size_t A, std::size_t... Bs>
  static constexpr bool one_against_all(std::index_sequence<Bs...>) {
    return (identity_sharing_ok<substructure_t, A, Bs>::value && ...);
  }

  template <std::size_t... As>
  static constexpr bool all_sharing_ok(std::index_sequence<As...>) {
    return (one_against_all<As>(std::make_index_sequence<cells>{}) && ...);
  }

  static constexpr bool value =
      all_constants_ok(std::make_index_sequence<cells>{}) &&
      all_sharing_ok(std::make_index_sequence<cells>{});
};

template <typename substructure_t>
struct is_unit_upper_triangular_helper<substructure_t, true> {
  template <std::size_t... Ks>
  static constexpr bool all_unit_diagonal(std::index_sequence<Ks...>) {
    constexpr std::size_t C = substructure_t::columns;
    return ((Ks / C != Ks % C ||
             std::same_as<
                 typename substructure_t::template value_at<Ks / C, Ks % C>,
                 one>) &&
            ...);
  }
  static constexpr bool value =
      is_upper_triangular_helper<substructure_t>::value &&
      all_unit_diagonal(std::make_index_sequence<substructure_t::rows *
                                                 substructure_t::columns>{});
};
}  // namespace detail

namespace concepts {
/// @brief True iff the type is a KSM row.
template <typename T>
concept is_row = detail::is_row_v<T>;

/// @brief True iff the type is a KSM row with only canonical values.
template <typename T>
concept is_canonical_row = detail::is_canonical_row_v<T>;

/// @brief True iff the type is a KSM substructure with only canonical values.
template <typename T>
concept is_canonical = detail::is_canonical_substructure_v<T>;

/// @brief True iff the type is a symmetric KSM substructure.
template <typename substructure_t>
concept is_symmetric = is_canonical<substructure_t> &&
                       detail::is_symmetric_helper<substructure_t>::value;

/// @brief True iff every cell off the diagonal is a structural zero. A
/// diagonal substructure is also a symmetric one.
template <typename substructure_t>
concept is_diagonal = is_canonical<substructure_t> &&
                      detail::is_diagonal_helper<substructure_t>::value;

/// @brief True iff the type is an upper triangular substructure, i.e. all
/// elements in the lower triangle are zero.
template <typename substructure_t>
concept is_upper_triangular =
    is_canonical<substructure_t> &&
    detail::is_upper_triangular_helper<substructure_t>::value;

/// @brief True iff the type is an upper triangular substructure, i.e. all
/// elements in the lower triangle are zero, and the diagonal is all ones.
template <typename substructure_t>
concept is_unit_upper_triangular =
    is_canonical<substructure_t> &&
    detail::is_unit_upper_triangular_helper<substructure_t>::value;

/// @brief True iff the substructure can hold the identity matrix.
///
/// Two things can stop it. A constant cell is fixed, so it has to already be
/// what the identity needs there; and cells sharing a variable share one
/// scalar, so they cannot disagree about it. A substructure naming the same
/// variable on and off the diagonal fails the second.
///
/// This does not require a square substructure: a rectangular one gets ones
/// along the leading diagonal and zeros elsewhere.
template <typename substructure_t>
concept can_represent_identity =
    is_canonical<substructure_t> &&
    detail::can_represent_identity_helper<substructure_t>::value;

/// @brief True iff multiplying a substructure by itself gives that same
/// substructure back.
///
/// A closed substructure is one a chain of products can stay inside: the type
/// of A*A is the type of A, so it is also the type of A*A*A. This is an
/// important class of substructures as it enables e.g. using the substructure
/// for propagating a Jacobian.
template <typename substructure_t>
concept is_closed_under_multiplication =
    is_canonical<substructure_t> &&
    (substructure_t::rows == substructure_t::columns) &&
    std::same_as<
        typename substructure_t::template multiplication_type<substructure_t>,
        substructure_t>;
}  // namespace concepts
}  // namespace detray::ksm
