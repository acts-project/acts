// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/algorithms/matrix/determinant/cofactor.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <cassert>
#include <type_traits>

namespace detray::algebra::generic::matrix {

namespace adjoint {

/// "Adjoint getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct cofactor {
  /// Index type used to traverse matrix rows and columns.
  using index_type = detray::traits::index_t<matrix_t>;

  /// Return the adjugate of a square matrix.
  /// @param m Matrix whose adjugate is calculated.
  /// @returns The transposed cofactor matrix.
  DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
    return adjoint_getter_helper<detray::traits::max_rank<matrix_t>>()(m);
  }

  /// Helper that computes the adjugate for a fixed matrix rank.
  template <index_type N>
  struct adjoint_getter_helper;

  /// Adjugate specialization for a one-element matrix.
  template <index_type N>
    requires(N == 1)
  struct adjoint_getter_helper<N> {
    /// Return the one-element identity matrix, the adjugate in rank one.
    /// @returns The one-element identity matrix.
    DETRAY_HOST_DEVICE constexpr matrix_t operator()(
        const matrix_t & /*m*/) const {
      using element_getter_t = detray::traits::element_getter_t<matrix_t>;
      constexpr element_getter_t elem{};

      matrix_t ret;
      elem(ret, 0, 0) = 1;

      return ret;
    }
  };

  /// Adjugate specialization that forms and transposes the cofactor matrix.
  template <index_type N>
    requires(N != 1)
  struct adjoint_getter_helper<N> {
    /// Determinant implementation used to evaluate each minor.
    using determinant_getter = determinant::cofactor<matrix_t>;

    /// Form the adjugate by computing signed minors and transposing them.
    /// @param m Matrix whose adjugate is calculated.
    /// @returns The transposed cofactor matrix.
    DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
      using index_t = detray::traits::index_t<matrix_t>;
      using element_getter_t = detray::traits::element_getter_t<matrix_t>;
      constexpr element_getter_t elem{};

      // temp is used to store cofactors of m
      int sign = 1;

      matrix_t temp;  //< To store cofactors
      matrix_t adj;

      for (index_t i = 0; i < N; i++) {
        for (index_t j = 0; j < N; j++) {
          // Get cofactor of m[i][j]
          typename determinant_getter::template determinant_getter_helper<N>()
              .get_cofactor(m, temp, i, j);

          // sign of adj[j][i] positive if sum of row
          // and column indexes is even.
          sign = ((i + j) % 2 == 0) ? 1 : -1;

          // Interchanging rows and columns to get the
          // transpose of the cofactor matrix
          elem(adj, j, i) =
              sign *
              typename determinant_getter::template determinant_getter_helper<
                  N - 1>()(temp);
        }
      }

      return adj;
    }
  };
};

}  // namespace adjoint

namespace inverse {

/// "inverse getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct cofactor {
  /// Scalar type of matrix elements.
  using scalar_t = detray::traits::value_t<matrix_t>;
  /// Index type used to traverse matrix rows and columns.
  using index_t = detray::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;
  /// Cofactor-expansion determinant implementation.
  using determinant_getter = determinant::cofactor<matrix_t>;
  /// Adjugate implementation used in the inverse formula.
  using adjoint_getter = adjoint::cofactor<matrix_t>;

  /// Invert a nonsingular square matrix.
  /// @param m Matrix to invert.
  /// @returns The inverse matrix; a singular input triggers an assertion.
  DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
    constexpr element_getter_t elem{};
    constexpr index_t N{detray::traits::max_rank<matrix_t>};

    matrix_t ret;

    // Find determinant of A
    scalar_t det = determinant_getter()(m);

    // TODO: handle singular matrix error
    if (det == 0) {
      assert(false && "Singular matrix during inversion (cofactor)");
      return ret;
    }

    auto adj = adjoint_getter()(m);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (index_t i = 0; i < N; i++) {
      for (index_t j = 0; j < N; j++) {
        elem(ret, j, i) = elem(adj, j, i) / det;
      }
    }

    return ret;
  }
};

}  // namespace inverse

}  // namespace detray::algebra::generic::matrix
