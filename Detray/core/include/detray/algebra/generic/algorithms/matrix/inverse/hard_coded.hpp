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
#include "detray/algebra/generic/algorithms/matrix/determinant/hard_coded.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <type_traits>

namespace detray::algebra::generic::matrix::inverse {

/// "inverse getter", assuming a N X N matrix
template <concepts::square_matrix matrix_t>
struct hard_coded {
  using scalar_t = detray::traits::value_t<matrix_t>;
  using index_t = detray::traits::index_t<matrix_t>;

  /// Function (object) used for accessing a matrix element
  using element_getter_t = detray::traits::element_getter_t<matrix_t>;
  using determinant_getter = determinant::hard_coded<matrix_t>;

  // 2 X 2 matrix inverse
  template <typename M = matrix_t>
    requires(detray::traits::max_rank<M> == 2)
  DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
    constexpr element_getter_t elem{};

    const scalar_t det = determinant_getter()(m);

    matrix_t ret;
    elem(ret, 0, 0) = elem(m, 1, 1) / det;
    elem(ret, 0, 1) = -elem(m, 0, 1) / det;
    elem(ret, 1, 0) = -elem(m, 1, 0) / det;
    elem(ret, 1, 1) = elem(m, 0, 0) / det;

    return ret;
  }

  // 4 X 4 matrix inverse
  template <typename M = matrix_t>
    requires(detray::traits::max_rank<M> == 4)
  DETRAY_HOST_DEVICE constexpr matrix_t operator()(const matrix_t &m) const {
    constexpr element_getter_t elem{};

    matrix_t ret;
    elem(ret, 0, 0) = elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 1) -
                      elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 1) +
                      elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 2) -
                      elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 2) -
                      elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 3) +
                      elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 3);
    elem(ret, 0, 1) = elem(m, 0, 3) * elem(m, 2, 2) * elem(m, 3, 1) -
                      elem(m, 0, 2) * elem(m, 2, 3) * elem(m, 3, 1) -
                      elem(m, 0, 3) * elem(m, 2, 1) * elem(m, 3, 2) +
                      elem(m, 0, 1) * elem(m, 2, 3) * elem(m, 3, 2) +
                      elem(m, 0, 2) * elem(m, 2, 1) * elem(m, 3, 3) -
                      elem(m, 0, 1) * elem(m, 2, 2) * elem(m, 3, 3);
    elem(ret, 0, 2) = elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 3, 1) -
                      elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 3, 1) +
                      elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 3, 2) -
                      elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 3, 2) -
                      elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 3, 3) +
                      elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 3, 3);
    elem(ret, 0, 3) = elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 1) -
                      elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 1) -
                      elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 2) +
                      elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 2) +
                      elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 3) -
                      elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 3);
    elem(ret, 1, 0) = elem(m, 1, 3) * elem(m, 2, 2) * elem(m, 3, 0) -
                      elem(m, 1, 2) * elem(m, 2, 3) * elem(m, 3, 0) -
                      elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 2) +
                      elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 2) +
                      elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 3) -
                      elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 3);
    elem(ret, 1, 1) = elem(m, 0, 2) * elem(m, 2, 3) * elem(m, 3, 0) -
                      elem(m, 0, 3) * elem(m, 2, 2) * elem(m, 3, 0) +
                      elem(m, 0, 3) * elem(m, 2, 0) * elem(m, 3, 2) -
                      elem(m, 0, 0) * elem(m, 2, 3) * elem(m, 3, 2) -
                      elem(m, 0, 2) * elem(m, 2, 0) * elem(m, 3, 3) +
                      elem(m, 0, 0) * elem(m, 2, 2) * elem(m, 3, 3);
    elem(ret, 1, 2) = elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 3, 0) -
                      elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 3, 0) -
                      elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 3, 2) +
                      elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 3, 2) +
                      elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 3, 3) -
                      elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 3, 3);
    elem(ret, 1, 3) = elem(m, 0, 2) * elem(m, 1, 3) * elem(m, 2, 0) -
                      elem(m, 0, 3) * elem(m, 1, 2) * elem(m, 2, 0) +
                      elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 2) -
                      elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 2) -
                      elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 3) +
                      elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 3);
    elem(ret, 2, 0) = elem(m, 1, 1) * elem(m, 2, 3) * elem(m, 3, 0) -
                      elem(m, 1, 3) * elem(m, 2, 1) * elem(m, 3, 0) +
                      elem(m, 1, 3) * elem(m, 2, 0) * elem(m, 3, 1) -
                      elem(m, 1, 0) * elem(m, 2, 3) * elem(m, 3, 1) -
                      elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 3) +
                      elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 3);
    elem(ret, 2, 1) = elem(m, 0, 3) * elem(m, 2, 1) * elem(m, 3, 0) -
                      elem(m, 0, 1) * elem(m, 2, 3) * elem(m, 3, 0) -
                      elem(m, 0, 3) * elem(m, 2, 0) * elem(m, 3, 1) +
                      elem(m, 0, 0) * elem(m, 2, 3) * elem(m, 3, 1) +
                      elem(m, 0, 1) * elem(m, 2, 0) * elem(m, 3, 3) -
                      elem(m, 0, 0) * elem(m, 2, 1) * elem(m, 3, 3);
    elem(ret, 2, 2) = elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 3, 0) -
                      elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 3, 0) +
                      elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 3, 1) -
                      elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 3, 1) -
                      elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 3, 3) +
                      elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 3, 3);
    elem(ret, 2, 3) = elem(m, 0, 3) * elem(m, 1, 1) * elem(m, 2, 0) -
                      elem(m, 0, 1) * elem(m, 1, 3) * elem(m, 2, 0) -
                      elem(m, 0, 3) * elem(m, 1, 0) * elem(m, 2, 1) +
                      elem(m, 0, 0) * elem(m, 1, 3) * elem(m, 2, 1) +
                      elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 3) -
                      elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 3);
    elem(ret, 3, 0) = elem(m, 1, 2) * elem(m, 2, 1) * elem(m, 3, 0) -
                      elem(m, 1, 1) * elem(m, 2, 2) * elem(m, 3, 0) -
                      elem(m, 1, 2) * elem(m, 2, 0) * elem(m, 3, 1) +
                      elem(m, 1, 0) * elem(m, 2, 2) * elem(m, 3, 1) +
                      elem(m, 1, 1) * elem(m, 2, 0) * elem(m, 3, 2) -
                      elem(m, 1, 0) * elem(m, 2, 1) * elem(m, 3, 2);
    elem(ret, 3, 1) = elem(m, 0, 1) * elem(m, 2, 2) * elem(m, 3, 0) -
                      elem(m, 0, 2) * elem(m, 2, 1) * elem(m, 3, 0) +
                      elem(m, 0, 2) * elem(m, 2, 0) * elem(m, 3, 1) -
                      elem(m, 0, 0) * elem(m, 2, 2) * elem(m, 3, 1) -
                      elem(m, 0, 1) * elem(m, 2, 0) * elem(m, 3, 2) +
                      elem(m, 0, 0) * elem(m, 2, 1) * elem(m, 3, 2);
    elem(ret, 3, 2) = elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 3, 0) -
                      elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 3, 0) -
                      elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 3, 1) +
                      elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 3, 1) +
                      elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 3, 2) -
                      elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 3, 2);
    elem(ret, 3, 3) = elem(m, 0, 1) * elem(m, 1, 2) * elem(m, 2, 0) -
                      elem(m, 0, 2) * elem(m, 1, 1) * elem(m, 2, 0) +
                      elem(m, 0, 2) * elem(m, 1, 0) * elem(m, 2, 1) -
                      elem(m, 0, 0) * elem(m, 1, 2) * elem(m, 2, 1) -
                      elem(m, 0, 1) * elem(m, 1, 0) * elem(m, 2, 2) +
                      elem(m, 0, 0) * elem(m, 1, 1) * elem(m, 2, 2);

    const scalar_t idet = static_cast<scalar_t>(1.) / determinant_getter()(ret);
    for (unsigned int c = 0; c < 4; ++c) {
      for (unsigned int r = 0; r < 4; ++r) {
        elem(ret, c, r) *= idet;
      }
    }

    return ret;
  }
};

}  // namespace detray::algebra::generic::matrix::inverse
