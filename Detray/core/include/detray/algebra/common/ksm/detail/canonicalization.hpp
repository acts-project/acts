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
#include "detray/algebra/common/ksm/detail/symbolic.hpp"
#include "detray/algebra/common/ksm/make_substructure.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

// This file defines the process of canonicalization, which is a process in
// which we turn a matrix substructure into a canonical format. This allows us
// to easily reason about structural properties of matrices even after complex
// series of operations.
namespace detray::ksm::detail {
/// @brief Internal helper function for `first_equal_flat_index_v`.
///
/// Iterates over a list of flat indices Ks and uses a constexpr loop to find
/// the first row-major flat index in Sub that matches Val.
template <typename Sub, typename Val, std::size_t... Ks>
constexpr std::size_t first_equal_flat_index_helper(
    std::index_sequence<Ks...>) {
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

/// @brief Canonicalization helper that finds the flat index of the first
/// matrix element that is equal to the one given.
///
/// @tparam Sub The substructure to search in (haystack).
/// @tparam Val The value to look for (needle).
///
/// Given a matrix substructure of:
///
/// [[A, B, C],
///  [D, B, A],
///  [D, E, A]]
///
/// A search for B would return number 1, as the first occurrence in a row-major
/// order is (0, 1) or flat index 1. Searching for D returns 3, as the first
/// index is (1, 0) or 3.
///
/// This wraps first_equal_flat_index_helper by automatically providing the
/// index sequence.
template <typename Sub, typename Val>
inline constexpr std::size_t first_equal_flat_index_v =
    first_equal_flat_index_helper<Sub, Val>(
        std::make_index_sequence<Sub::rows * Sub::columns>{});

/// @brief Helper function for canonicalization, which determines the canonical
/// value of a matrix at a given coordinate.
///
/// @tparam T The substructure to be canonicalized.
/// @tparam I The row index to compute the canonical element for.
/// @tparam J The column index to compute the canonical element for.
///
/// See @c canonicalize_substructure for more details on how canonicalization
/// works.
template <typename T, std::size_t I, std::size_t J>
struct canonicalize_substructure_helper {
 private:
  using curr_value = T::template value_at<I, J>;
  using resolved_value = typename resolve_value<curr_value>::type;

 public:
  // Encoded here is the core canonicalization logic. For any given (I, J)
  // coordinate, we do the following:
  //
  // 1. Check if the value at (I, J) is an anonymous variable. If yes, return
  //    a non-anonymous value with flat index (I, J).
  // 2. Otherwise, check if the value is a constant integral value. If yes,
  //    this value is already canonical and we return it as-is.
  // 3. Otherwise, find the flat index of the first matrix element equal to
  //    the element at (I, J) (which might well be (I, J)) and return a
  //    non-anonymous variable with that index.
  using type = std::conditional_t<
      std::same_as<curr_value, variable>, index_variable<T::columns * I + J>,
      std::conditional_t<
          concepts::is_integral_value<resolved_value>, resolved_value,
          index_variable<first_equal_flat_index_v<T, curr_value>>>>;
};

/// @brief A type-level function canonicalize a matrix substructure.
///
/// A matrix substructure is canonical if:
///
/// 1. The substructure contains only non-anonymous variables or constant
///    values; and
/// 2. The application of a chosen canonicalization function is idempotent on
///    it.
///
/// Where a canonicalization function is a function that:
///
/// 1. Produces a canonical representation when passed any substructure.
/// 2. Maps structurally identical matrix values to the same value.
/// 3. Does not map structurally dissimilar values to the same value.
///
/// Note that there are very many possible canonicalization functions which
/// meet the aforementioned properties, this is just one. This is also why the
/// aforementioned definition of canonicality is circular.
///
/// The canonicalization function here works by mapping constant values to
/// constant values (which a canonicalization function is definitionally
/// allowed to do), and it maps any other value to a non-anonymous variable at
/// the index where the _first_ identical value is found. Note that for this
/// metric, an anonymous value is _never_ identical to any other value.
///
/// Take, as an example, the following matrix (where 0..9 are constants, _ is
/// an anonymous variable, and A..Z are non-anonymous variables). This matrix
/// substructure is non-canonical.
///
/// [[   A+1,     A,     0],
///  [   A+B,     C,     _],
///  [     1,     B,   A+B]]
///
/// In the canonical representation, we will refer to non-anonymous variables
/// using a tuple IJ=X, with the I, J coordinates and the corresponding flat
/// coordinate.
///
/// The canonical substructure for the matrix above is:
///
/// [[ 0,0=0, 0,1=1,     0],
///  [ 1,0=3, 1,1=4, 1,2=5],
///  [     1, 2,1=7, 1,0=3]]
///
/// Translating the non-anonymous indices back to letter forms as [0,1,2,...] =
/// [A,B,C,...] gives:
///
/// [[     A,     B,     0],
///  [     D,     E,     F],
///  [     1,     H,     D]]
///
/// Note that this procedure has preserved constants, and has mapped both cells
/// with value A+B to a new non-anonymous variable D.
///
/// Confirming that the properties of canonicalization hold is left as an
/// exercise to the reader. ;)
template <typename T>
struct canonicalize_substructure {
  // We curry the canonicalize_substructure_helper function with the matrix
  // substructure so that we can call generate_substructure_rows on it.
  template <std::size_t I, std::size_t J>
  using accessor = canonicalize_substructure_helper<T, I, J>;

  // To compute the canonical substructure, simply map over a matrix of the
  // same shape and call into the canonicalize substructure helper.
  using type = typename generate_substructure_rows<
      accessor, T::columns, std::make_index_sequence<T::rows>>::type;
};
}  // namespace detray::ksm::detail
