// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_aos_concepts.hpp"
#include "algebra/impl/vc_aos_vector.hpp"
#include "detray/algebra/common/impl/matrix.hpp"
#include "detray/algebra/common/storage/matrix.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/generic.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::vc_aos::math {

using algebra::common::identity;
using algebra::common::set_identity;
using algebra::common::set_zero;
using algebra::common::transpose;
using algebra::common::zero;

/// @returns the determinant
template <std::size_t N, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr value_t determinant(
    const algebra::storage::matrix<array_t, value_t, N, N> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::determinant(m);
}

/// @returns the inverse
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, value_t, COL,
                                                      ROW>
inverse(
    const algebra::storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::inverse(m);
}

/// @returns the transpose
template <std::size_t ROW, std::size_t COL, concepts::value value_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, value_t, COL,
                                                      ROW>
transpose(
    const algebra::storage::matrix<array_t, value_t, ROW, COL> &m) noexcept {
  // TODO: Implement vectorization friendly version
  return algebra::generic::math::transpose(m);
}

/// Column-wise cross product between a 3x3 matrix and a 3D vector
///
/// @note Uses the Vc AoS cross product, so that the result is the same as for
/// @c vector::cross on the matrix columns
///
/// @tparam value_t Value type of the matrix elements
/// @tparam array_t Underlying storage array template
/// @tparam vector_t Vc AoS 3D vector type
///
/// @param m the input matrix
/// @param v the input vector
///
/// @returns the matrix with the cross product of every column of @p m with
/// @p v
template <concepts::value value_t,
          template <typename, std::size_t> class array_t,
          concepts::vc_aos_vector vector_t>
  requires concepts::vector3D<vector_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::matrix<array_t, value_t, 3, 3>
column_wise_cross(const algebra::storage::matrix<array_t, value_t, 3, 3> &m,
                  const vector_t &v) {
  using matrix_t = algebra::storage::matrix<array_t, value_t, 3, 3>;
  using block_getter_t = detray::traits::block_getter_t<matrix_t>;
  constexpr block_getter_t block{};

  matrix_t ret{};

  for (std::size_t j = 0u; j < 3u; ++j) {
    const vector_t m_col = block.template vector<3>(m, 0u, j);
    block.set(ret,
              static_cast<vector_t>(algebra::vc_aos::math::cross(m_col, v)), 0u,
              j);
  }

  return ret;
}

}  // namespace detray::algebra::vc_aos::math
