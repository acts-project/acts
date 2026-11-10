// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/detail/math_helpers.hpp"
#include "detray/algebra/common/ksm/detail/symbolic.hpp"
#include "detray/algebra/common/ksm/make_substructure.hpp"
#include "detray/algebra/common/ksm/row.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <cstddef>
#include <type_traits>
#include <utility>

// This file encodes the substructure of matrices that are the result of
// certain operations like addition, subtraction, etc.
//
// All of these functions only work on canonical substructures. Always
// canonicalize first!
namespace detray::ksm::detail {
/// @brief Type function that computes element (I, J) of the addition of
/// substructures A+B.
template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(concepts::is_canonical<substructure_A> &&
           concepts::is_canonical<substructure_B>)
struct get_addition_element {
  using A_value_type = substructure_A::template value_at<I, J>;
  using B_value_type = substructure_B::template value_at<I, J>;

  // The result is simply defined as A_IJ + B_IJ.
  using type = value_addition<A_value_type, B_value_type>;
};

/// @brief Type function that computes element (I, J) of the subtraction of
/// substructures A-B.
template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(concepts::is_canonical<substructure_A> &&
           concepts::is_canonical<substructure_B>)
struct get_subtraction_element {
  using A_value_type = substructure_A::template value_at<I, J>;
  using B_value_type = substructure_B::template value_at<I, J>;

  // The result is simply defined as A_IJ - B_IJ.
  using type = value_subtraction<A_value_type, B_value_type>;
};

/// @brief Type functions to compute a dot product of two rows.
///
/// Used to define matrix multiplication.
///
/// @{
template <typename T, typename S>
struct dot_product {};

// The obvious first specialization here is for non-empty rows.
template <typename T, typename... Ts, typename S, typename... Ss>
struct dot_product<row<T, Ts...>, row<S, Ss...>> {
  // Dot product is defined only for rows of equal size.
  static_assert(sizeof...(Ts) == sizeof...(Ss));

 private:
  // We compute both the value at the head of the rows (i.e. the product of
  // the first value in each row) and then we recurse into this function to
  // get the remainder.
  //
  // NOTE: This computation is technically less strict than it could be, as
  // it collapses the current value and remaining values to variables.
  // However, there is relatively little structural advantage to be gained in
  // matrix multiplication.
  using rem_val = typename dot_product<row<Ts...>, row<Ss...>>::type;
  using cur_val = decltype(resolve_symbolic<value_multiplication<T, S>>());

 public:
  // The resulting type is simply the resolved symbolic computation of the
  // addition of the head and tail values.
  using type = decltype(resolve_symbolic<value_addition<cur_val, rem_val>>());
};

// Then, as the dot product is the sum, the base case is simply the additive
// identity, zero.
template <>
struct dot_product<row<>, row<>> {
  using type = zero;
};
/// @}

/// @brief Type function that computes element (I, J) of the multiplication of
/// substructures A*B.
///
/// This is somewhat more complicated than addition or subtraction, as it needs
/// a dot product!
template <typename substructure_A, typename substructure_B, std::size_t I,
          std::size_t J>
  requires(concepts::is_canonical<substructure_A> &&
           concepts::is_canonical<substructure_B>)
struct get_multiplication_element {
 private:
  using A_row_type = typename substructure_A::template row_at<I>;
  using B_column_type = typename substructure_B::template column_at<J>;

  static_assert(A_row_type::num_values == B_column_type::num_values);

 public:
  using type = typename dot_product<A_row_type, B_column_type>::type;
};

/// @brief Type function that computes element (I, J) of a matrix A after
/// applying the identity to it, i.e. just the element A_IJ!
template <typename substructure_A, std::size_t I, std::size_t J>
  requires(concepts::is_canonical<substructure_A>)
struct get_identity_element {
  using type = typename substructure_A::template value_at<I, J>;
};

/// @brief Type function that computes element (I, J) of a matrix A after
/// negating it, i.e. -A_IJ.
template <typename substructure_A, std::size_t I, std::size_t J>
  requires(concepts::is_canonical<substructure_A>)
struct get_negation_element {
  using type = value_negation<typename substructure_A::template value_at<I, J>>;
};

/// @brief Type function that computes the substructure of the addition of
/// substructures A+B.
template <typename A, typename B>
  requires(concepts::is_canonical<A> && concepts::is_canonical<B>)
struct additive_substructure {
  static_assert(A::columns == B::columns);
  static_assert(A::rows == B::rows);

  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  // We simply delegate all logic into a curried version of
  // get_addition_element!
  template <std::size_t I, std::size_t J>
  using accessor = get_addition_element<A, B, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

/// @brief Type function that computes the substructure of the subtraction of
/// substructures A-B.
template <typename A, typename B>
  requires(concepts::is_canonical<A> && concepts::is_canonical<B>)
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

/// @brief Type function that computes the substructure of the multiplication of
/// substructures A*B.
template <typename A, typename B>
  requires(concepts::is_canonical<A> && concepts::is_canonical<B>)
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

/// @brief Type function that computes the substructure of negated substructure
/// -A.
template <typename A>
  requires(concepts::is_canonical<A>)
struct negated_substructure {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  template <std::size_t I, std::size_t J>
  using accessor = get_negation_element<A, I, J>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

/// @brief Type function that computes the substructure of A, but with the
/// upper triangle copied into the lower triangle.
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

/// @brief Helper function that conditionally upper-triangularizes a
/// substructure.
///
/// Helpful for logic like "produce a symmetric substructure only if the input
/// to a function is symmetric".
///
/// @{
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
/// @}

/// @brief Type functions that compute the substructure of inverted
/// substructure A^-1.
///
/// @note This is only supported for dimensions 1 and 2 for now.
template <typename A>
struct inverted_substructure {};

// The 1D case.
template <typename A>
  requires(concepts::is_canonical<A> && concepts::is_invertible<A> &&
           A::rows == 1)
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

  // Note that symmetrizing a 1D matrix makes little sense, although it also
  // does no harm...
  using type = upper_symmetric_if<raw_type, concepts::is_symmetric<A>>::type;
};

// The 2D case.
template <typename A>
  requires(concepts::is_canonical<A> && concepts::is_invertible<A> &&
           A::rows == 2)
struct inverted_substructure<A> {
  static constexpr std::size_t columns = A::columns;
  static constexpr std::size_t rows = A::rows;

  using det = get_determinant_helper<A>::type;

  template <std::size_t I, std::size_t J>
  struct accessor {
    // This just models the standard formula for 2x2 inversion.
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

  // Here, we explicitly symmetrize the inverse of a symmetric 2D matrix, as
  // this symmetry is difficult (but not impossible) to extract through
  // canonicalization.
  using type = upper_symmetric_if<raw_type, concepts::is_symmetric<A>>::type;
};

/// @brief Type functions that compute the substructure of transposed
/// substructure A^T.
template <typename A>
  requires(concepts::is_canonical<A>)
struct transposed_substructure {
  static constexpr std::size_t columns = A::rows;
  static constexpr std::size_t rows = A::columns;

  // Here the identity element comes in handy, we just swap the arguments!
  //
  // That's functional programming for you.
  template <std::size_t I, std::size_t J>
  using accessor = get_identity_element<A, J, I>;

  using type =
      typename generate_substructure_rows<accessor, columns,
                                          std::make_index_sequence<rows>>::type;
};

/// @brief Type functions that compute the substructure of congruence
/// substructure ABA^T.
template <typename A, typename B>
  requires(concepts::is_canonical<A> && concepts::is_canonical<B>)
struct congruence_substructure {
  static_assert(A::columns == B::rows);
  static_assert(B::columns == A::columns);

  using AB_type =
      typename multiplicative_substructure<A, B>::type::canonical_type;
  using ABAT_type = typename multiplicative_substructure<
      AB_type, typename transposed_substructure<A>::type::canonical_type>::
      type::canonical_type;

  // This substructure is also explicitly symmetric iff. B is symmetric.
  //
  // Proof:
  //
  // Q is symmetric if Q = Q^T
  //
  // ABA^T = (ABA^T)^T
  //       = ((A^T)^T)(B^T)(A^T)
  //       = AB^TA^T
  //
  // Iff B is symmetric, B = B^T so:
  //
  // AB^TA^T = ABA^T
  using sym_ABAT_type =
      typename upper_symmetric_substructure<ABAT_type>::type::canonical_type;
  static_assert(concepts::is_symmetric<sym_ABAT_type>);

  static constexpr std::size_t columns = A::rows;
  static constexpr std::size_t rows = A::rows;

  using type = sym_ABAT_type;
};
}  // namespace detray::ksm::detail
