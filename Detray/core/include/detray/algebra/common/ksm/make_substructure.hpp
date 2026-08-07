// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/fwd.hpp"
#include "detray/algebra/common/ksm/value.hpp"

// System include(s)
#include <cstddef>
#include <type_traits>
#include <utility>

// This file contains logic to create some common substructures.
namespace detray::ksm {
namespace detail {
/// @brief The bread and butter of all substructure generation, these types
/// respectively generate a row and a whole substructure using a callable.
///
/// @{
///
/// @brief Helper functions for @c generate_substructure_rows
///
/// @{
template <template <std::size_t, std::size_t> typename C, std::size_t I,
          typename T>
struct generate_substructure_values {};

template <template <std::size_t, std::size_t> typename C, std::size_t I,
          std::size_t... Js>
struct generate_substructure_values<C, I,
                                    std::integer_sequence<std::size_t, Js...>> {
  using type = row<typename C<I, Js>::type...>;
};
/// @}

/// @brief Generate a substructure using a generator function and a size.
///
/// This type function takes a function (I, J) -> Value which generates the
/// value at each index, then it takes a number of columns N, and a integer
/// sequence of row indices Is...
///
/// @{
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
/// @}
/// @}

/// @brief Helper function for @c symmetric_cell; computes the packed index of
/// a cell of an NxN symmetric matrix, counting only the upper triangle.
///
/// A cell in the lower triangle is mirrored into the upper one first, so (i, j)
/// and (j, i) get the same index, which is what makes them share a variable.
/// The indices run 0 to N*(N+1)/2 - 1 in row-major order *within the upper
/// triangle*, so they are not row-major indices of the full matrix: for N = 3,
/// (1, 1) is 3 here and 4 there.
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

/// @brief Helper for fully dense matrices; always returns an anonymous
/// variable.
template <std::size_t I, std::size_t J>
struct dense_cell {
  using type = variable;  // every cell an independent variable
};

/// @brief Helper for diagonal matrices; variable on the diagonal and zero
/// elsewhere.
template <std::size_t I, std::size_t J>
struct diagonal_cell {
  using type = std::conditional_t<(I == J), variable, zero>;
};

/// @brief Helper for identity matrices; one on the diagonal and zero
/// elsewhere.
template <std::size_t I, std::size_t J>
struct identity_cell {
  using type = std::conditional_t<(I == J), one, zero>;
};

/// @brief Helper for symmetric matrices. A cell in the lower triangle names
/// the same variable as its mirror in the upper triangle, so the two share one
/// stored scalar.
template <std::size_t N>
struct symmetric_cell {
  template <std::size_t I, std::size_t J>
  struct cell {
    using type = index_variable<symmetric_index<N>(I, J)>;
  };
};

/// @brief Helper for matrices whose K leftmost columns hold variables, every
/// further column being a structural zero.
template <std::size_t K>
struct left_columns_cell {
  template <std::size_t I, std::size_t J>
  struct cell {
    using type = std::conditional_t<(J < K), variable, zero>;
  };
};

}  // namespace detail

/// @brief General substructure creator. Takes a helper function, a number
/// of rows and columns, and generates a substructure.
template <template <std::size_t, std::size_t> typename C, std::size_t Rows,
          std::size_t Cols>
using make_substructure = typename detail::generate_substructure_rows<
    C, Cols, std::make_integer_sequence<std::size_t, Rows>>::type;

/// @brief Generate a fully dense substructure.
template <std::size_t Rows, std::size_t Cols = Rows>
using make_dense_substructure =
    make_substructure<detail::dense_cell, Rows, Cols>;

/// @brief Generate a symmetric substructure.
template <std::size_t N>
using make_symmetric_substructure =
    make_substructure<detail::symmetric_cell<N>::template cell, N, N>;

/// @brief Generate a diagonal substructure.
template <std::size_t N>
using make_diagonal_substructure =
    make_substructure<detail::diagonal_cell, N, N>;

/// @brief Generate an identity substructure.
template <std::size_t N>
using make_identity_substructure =
    make_substructure<detail::identity_cell, N, N>;

/// @brief Generate a substructure whose K leftmost columns hold variables,
/// every further column being a structural zero.
template <std::size_t Rows, std::size_t Cols, std::size_t K>
using make_left_columns_substructure =
    make_substructure<detail::left_columns_cell<K>::template cell, Rows, Cols>;
}  // namespace detray::ksm
