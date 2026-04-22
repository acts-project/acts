// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace detray::algebra::generic::matrix::determinant {

/// "Determinant getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct cofactor {
  using scalar_type = detray::traits::value_t<matrix_t>;
  using index_type = detray::traits::index_t<matrix_t>;

  DETRAY_HOST_DEVICE constexpr scalar_type operator()(const matrix_t &m) const {
    return determinant_getter_helper<detray::traits::max_rank<matrix_t>>()(m);
  }

  template <index_type N>
  struct determinant_getter_helper;

  template <index_type N>
    requires(N == 1)
  struct determinant_getter_helper<N> {
    template <class input_matrix_type>
    DETRAY_HOST_DEVICE constexpr scalar_type operator()(
        const input_matrix_type &m) const {
      using element_getter_t =
          detray::traits::element_getter_t<input_matrix_type>;

      constexpr element_getter_t elem{};

      return elem(m, 0, 0);
    }
  };

  template <index_type N>
    requires(N != 1)
  struct determinant_getter_helper<N> {
    template <class input_matrix_type>
    DETRAY_HOST_DEVICE constexpr scalar_type operator()(
        const input_matrix_type &m) const {
      using scalar_t = detray::traits::value_t<input_matrix_type>;
      using index_t = detray::traits::index_t<input_matrix_type>;
      using element_getter_t =
          detray::traits::element_getter_t<input_matrix_type>;

      constexpr element_getter_t elem{};

      scalar_t D = 0;

      // To store cofactors
      input_matrix_type temp;

      // To store sign multiplier
      int sign = 1;

      // Iterate for each element of first row
      for (index_t col = 0; col < N; col++) {
        // Getting Cofactor of A[0][f]
        this->get_cofactor(m, temp, index_t(0), col);
        D += sign * elem(m, 0, col) * determinant_getter_helper<N - 1>()(temp);

        // terms are to be added with alternate sign
        sign = -sign;
      }

      return D;
    }

    template <class input_matrix_type>
    DETRAY_HOST_DEVICE constexpr void get_cofactor(const input_matrix_type &m,
                                                   input_matrix_type &temp,
                                                   index_type p,
                                                   index_type q) const {
      using index_t = detray::traits::index_t<input_matrix_type>;
      using element_getter_t =
          detray::traits::element_getter_t<input_matrix_type>;

      constexpr element_getter_t elem{};

      index_t i = 0;
      index_t j = 0;

      // Looping for each element of the matrix
      for (index_t row = 0; row < N; row++) {
        for (index_t col = 0; col < N; col++) {
          //  Copying into temporary matrix only those element
          //  which are not in given row and column
          if (row != p && col != q) {
            elem(temp, i, j++) = elem(m, row, col);

            // Row is filled, so increase row index and
            // reset col index
            if (j == N - 1) {
              j = 0;
              i++;
            }
          }
        }
      }
    }
  };
};

}  // namespace detray::algebra::generic::matrix::determinant
