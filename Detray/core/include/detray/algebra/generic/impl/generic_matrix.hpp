// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/algorithms/utils/algorithm_selector.hpp"
#include "detray/algebra/generic/impl/generic_vector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::generic::math {

/// Create zero matrix - generic transform3
template <concepts::matrix M>
DETRAY_HOST_DEVICE constexpr M zero() {
  using index_t = detray::traits::index_t<M>;
  using element_getter_t = detray::traits::element_getter_t<M>;

  M ret;

  for (index_t j = 0; j < detray::traits::columns<M>; ++j) {
    for (index_t i = 0; i < detray::traits::rows<M>; ++i) {
      element_getter_t{}(ret, i, j) = 0;
    }
  }

  return ret;
}

/// Create identity matrix - generic transform3
template <concepts::matrix M>
DETRAY_HOST_DEVICE constexpr M identity() {
  using index_t = detray::traits::index_t<M>;
  using element_getter_t = detray::traits::element_getter_t<M>;

  auto ret{zero<M>()};

  for (index_t i = 0; i < detray::traits::max_rank<M>; ++i) {
    element_getter_t{}(ret, i, i) = 1;
  }

  return ret;
}

/// Set @param m as zero matrix
template <concepts::matrix M>
DETRAY_HOST_DEVICE constexpr void set_zero(M &m) {
  m = zero<M>();
}

/// Set @param m as identity matrix
template <concepts::matrix M>
DETRAY_HOST_DEVICE constexpr void set_identity(M &m) {
  m = identity<M>();
}

/// @returns the transpose matrix of @param m
template <concepts::matrix M>
DETRAY_HOST_DEVICE constexpr auto transpose(const M &m) {
  using index_t = detray::traits::index_t<M>;
  using value_t = detray::traits::value_t<M>;
  using element_getter_t = detray::traits::element_getter_t<M>;

  constexpr index_t rows{detray::traits::rows<M>};
  constexpr index_t columns{detray::traits::columns<M>};

  detray::traits::get_matrix_t<M, columns, rows, value_t> ret;

  for (index_t i = 0; i < rows; ++i) {
    for (index_t j = 0; j < columns; ++j) {
      element_getter_t{}(ret, j, i) = element_getter_t{}(m, i, j);
    }
  }

  return ret;
}

/// Column-wise cross product between matrix (m) and vector (v)
template <concepts::square_matrix M, concepts::vector3D V>
  requires(detray::traits::max_rank<M> == 3u)
DETRAY_HOST_DEVICE constexpr M column_wise_cross(const M &m, const V &v) {
  using block_getter_t = detray::traits::block_getter_t<M>;
  constexpr block_getter_t block{};

  M ret;

  const V m_col0 = block.template vector<3>(m, 0u, 0u);
  const V m_col1 = block.template vector<3>(m, 0u, 1u);
  const V m_col2 = block.template vector<3>(m, 0u, 2u);

  block.set(ret, static_cast<V>(cross(m_col0, v)), 0u, 0u);
  block.set(ret, static_cast<V>(cross(m_col1, v)), 0u, 1u);
  block.set(ret, static_cast<V>(cross(m_col2, v)), 0u, 2u);

  return ret;
}

/// Column-wise multiplication between matrix (m) and vector (v)
template <concepts::square_matrix M, concepts::vector V>
  requires(detray::traits::max_rank<M> == detray::traits::size<V>)
DETRAY_HOST_DEVICE inline M column_wise_multiply(const M &m, const V &v) {
  using element_getter_t = detray::traits::element_getter_t<M>;
  constexpr element_getter_t elem{};

  M ret;

  constexpr std::size_t N{detray::traits::max_rank<M>};
  for (std::size_t i = 0u; i < N; i++) {
    for (std::size_t j = 0u; j < N; j++) {
      elem(ret, j, i) = elem(m, j, i) * elem(v, j);
    }
  }

  return ret;
}

/// Cross product matrix
template <concepts::vector3D V3>
DETRAY_HOST_DEVICE inline auto cross_matrix(const V3 &v) {
  using scalar_t = detray::traits::scalar_t<V3>;
  using matrix_33_t = detray::traits::get_matrix_t<V3, 3, 3, scalar_t>;

  using element_getter_t = detray::traits::element_getter_t<V3>;
  constexpr element_getter_t elem{};

  matrix_33_t ret;
  elem(ret, 0u, 0u) = 0.f;
  elem(ret, 0u, 1u) = -elem(v, 2u);
  elem(ret, 0u, 2u) = elem(v, 1u);
  elem(ret, 1u, 0u) = elem(v, 2u);
  elem(ret, 1u, 1u) = 0.f;
  elem(ret, 1u, 2u) = -elem(v, 0u);
  elem(ret, 2u, 0u) = -elem(v, 1u);
  elem(ret, 2u, 1u) = elem(v, 0u);
  elem(ret, 2u, 2u) = 0.f;

  return ret;
}

/// Outer product operation
template <concepts::vector V>
DETRAY_HOST_DEVICE inline detray::traits::get_matrix_t<
    V, detray::traits::size<V>, detray::traits::size<V>,
    detray::traits::scalar_t<V>>
outer_product(const V &v1, const V &v2) {
  using scalar_t = detray::traits::scalar_t<V>;
  using index_t = detray::traits::index_t<V>;

  constexpr auto N{static_cast<index_t>(detray::traits::size<V>)};
  using matrix_t = detray::traits::get_matrix_t<V, N, N, scalar_t>;

  using element_getter_t = detray::traits::element_getter_t<V>;
  constexpr element_getter_t elem{};

  matrix_t ret;

  // TODO: vectorize better
  for (index_t row = 0; row < N; row++) {
    for (index_t col = 0; col < N; col++) {
      elem(ret, row, col) = elem(v1, row) * elem(v2, col);
    }
  }

  return ret;
}

// Set matrix C to the product AB
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_product(MC &C, const MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable_into<MA, MB, MC>)
{
  using index_t = detray::traits::index_t<MC>;
  using value_t = detray::traits::value_t<MC>;

  for (index_t i = 0; i < detray::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < detray::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < detray::traits::rows<MB>; ++k) {
        t += detray::traits::element_getter_t<MA>()(A, i, k) *
             detray::traits::element_getter_t<MB>()(B, k, j);
      }

      detray::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix C to the product A^TB
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_product_left_transpose(MC &C, const MA &A,
                                                             const MB &B)
  requires(detray::concepts::matrix_multipliable_into<
           decltype(transpose(std::declval<MA>())), MB, MC>)
{
  using index_t = detray::traits::index_t<MC>;
  using value_t = detray::traits::value_t<MC>;

  for (index_t i = 0; i < detray::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < detray::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < detray::traits::rows<MB>; ++k) {
        t += detray::traits::element_getter_t<MA>()(A, k, i) *
             detray::traits::element_getter_t<MB>()(B, k, j);
      }

      detray::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix C to the product AB^T
template <concepts::matrix MC, concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_product_right_transpose(MC &C,
                                                              const MA &A,
                                                              const MB &B)
  requires(detray::concepts::matrix_multipliable_into<
           MA, decltype(transpose(std::declval<MB>())), MC>)
{
  using index_t = detray::traits::index_t<MC>;
  using value_t = detray::traits::value_t<MC>;

  for (index_t i = 0; i < detray::traits::rows<MC>; ++i) {
    for (index_t j = 0; j < detray::traits::columns<MC>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < detray::traits::columns<MA>; ++k) {
        t += detray::traits::element_getter_t<MA>()(A, i, k) *
             detray::traits::element_getter_t<MB>()(B, j, k);
      }

      detray::traits::element_getter_t<MC>()(C, i, j) = t;
    }
  }
}

// Set matrix A to the product AB in place
template <concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_inplace_product_right(MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable_into<MA, MB, MA>)
{
  using index_t = detray::traits::index_t<MA>;
  using value_t = detray::traits::value_t<MA>;

  for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
    detray::traits::get_matrix_t<MA, 1, detray::traits::columns<MA>, value_t> Q;

    for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
      detray::traits::element_getter_t<decltype(Q)>()(Q, 0, j) =
          detray::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
      value_t t = 0.f;

      for (index_t k = 0; k < detray::traits::rows<MB>; ++k) {
        t += detray::traits::element_getter_t<decltype(Q)>()(Q, 0, k) *
             detray::traits::element_getter_t<MB>()(B, k, j);
      }

      detray::traits::element_getter_t<MA>()(A, i, j) = t;
    }
  }
}

// Set matrix A to the product BA in place
template <concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_inplace_product_left(MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable_into<MB, MA, MA>)
{
  using index_t = detray::traits::index_t<MA>;
  using value_t = detray::traits::value_t<MA>;

  for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
    detray::traits::get_matrix_t<MA, 1, detray::traits::columns<MA>, value_t> Q;

    for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
      detray::traits::element_getter_t<decltype(Q)>()(Q, 0, i) =
          detray::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
      value_t t = 0.f;

      for (index_t k = 0; k < detray::traits::columns<MB>; ++k) {
        t += detray::traits::element_getter_t<MB>()(B, i, k) *
             detray::traits::element_getter_t<decltype(Q)>()(Q, 0, k);
      }

      detray::traits::element_getter_t<MA>()(A, i, j) = t;
    }
  }
}

// Set matrix A to the product AB^T in place
template <concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_inplace_product_right_transpose(
    MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable_into<
           MA, decltype(transpose(std::declval<MB>())), MA>)
{
  using index_t = detray::traits::index_t<MA>;
  using value_t = detray::traits::value_t<MA>;

  for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
    detray::traits::get_matrix_t<MA, 1, detray::traits::columns<MA>, value_t> Q;

    for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
      detray::traits::element_getter_t<decltype(Q)>()(Q, 0, j) =
          detray::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
      value_t T = 0.f;

      for (index_t k = 0; k < detray::traits::columns<MB>; ++k) {
        T += detray::traits::element_getter_t<decltype(Q)>()(Q, 0, k) *
             detray::traits::element_getter_t<MB>()(B, j, k);
      }

      detray::traits::element_getter_t<MA>()(A, i, j) = T;
    }
  }
}

// Set matrix A to the product B^TA in place
template <concepts::matrix MA, concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr void set_inplace_product_left_transpose(
    MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable_into<
           decltype(transpose(std::declval<MB>())), MA, MA>)
{
  using index_t = detray::traits::index_t<MA>;
  using value_t = detray::traits::value_t<MA>;

  for (index_t j = 0; j < detray::traits::columns<MA>; ++j) {
    detray::traits::get_matrix_t<MA, 1, detray::traits::columns<MA>, value_t> Q;

    for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
      detray::traits::element_getter_t<decltype(Q)>()(Q, 0, i) =
          detray::traits::element_getter_t<MA>()(A, i, j);
    }

    for (index_t i = 0; i < detray::traits::rows<MA>; ++i) {
      value_t T = 0.f;

      for (index_t k = 0; k < detray::traits::rows<MB>; ++k) {
        T += detray::traits::element_getter_t<MB>()(B, k, i) *
             detray::traits::element_getter_t<decltype(Q)>()(Q, 0, k);
      }

      detray::traits::element_getter_t<MA>()(A, i, j) = T;
    }
  }
}

template <bool transpose, concepts::matrix M>
DETRAY_HOST_DEVICE constexpr detray::traits::value_t<M> transposable_get(
    const M &A, const detray::traits::index_t<M> i,
    const detray::traits::index_t<M> j) {
  if constexpr (transpose) {
    return detray::traits::element_getter_t<M>()(A, j, i);
  } else {
    return detray::traits::element_getter_t<M>()(A, i, j);
  }
}

template <bool transpose_left, bool transpose_right, concepts::matrix MA,
          concepts::matrix MB>
DETRAY_HOST_DEVICE constexpr auto transposed_product(const MA &A, const MB &B)
  requires(detray::concepts::matrix_multipliable<
           std::conditional_t<transpose_left,
                              decltype(transpose(std::declval<MA>())), MA>,
           std::conditional_t<transpose_right,
                              decltype(transpose(std::declval<MB>())), MB>>)
{
  using index_t = detray::traits::index_t<MA>;
  using value_t = detray::traits::value_t<MA>;

  using effective_lhs_t =
      std::conditional_t<transpose_left,
                         decltype(transpose(std::declval<MA>())), MA>;
  using effective_rhs_t =
      std::conditional_t<transpose_right,
                         decltype(transpose(std::declval<MB>())), MB>;

  using result_t =
      detray::traits::get_matrix_t<MA, detray::traits::rows<effective_lhs_t>,
                                   detray::traits::columns<effective_rhs_t>,
                                   value_t>;

  result_t C = zero<result_t>();

  for (index_t i = 0; i < detray::traits::columns<effective_lhs_t>; ++i) {
    for (index_t j = 0; j < detray::traits::columns<result_t>; ++j) {
      for (index_t k = 0; k < detray::traits::rows<result_t>; ++k) {
        detray::traits::element_getter_t<result_t>()(C, k, j) +=
            transposable_get<transpose_left>(A, k, i) *
            transposable_get<transpose_right>(B, i, j);
      }
    }
  }

  return C;
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<M> determinant(
    const M &m) {
  return determinant_t<M>{}(m);
}

/// @returns the determinant of @param m
template <concepts::square_matrix M>
DETRAY_HOST_DEVICE constexpr M inverse(const M &m) {
  return inversion_t<M>{}(m);
}

/// @returns the fatcor L in the decomposition of @param m
template <concepts::square_matrix M>
DETRAY_HOST_DEVICE constexpr M cholesky_decomposition(const M &m) {
  using index_t = detray::traits::index_t<M>;
  using scalar_t = detray::traits::scalar_t<M>;
  using element_getter_t = detray::traits::element_getter_t<M>;

  constexpr element_getter_t elem{};
  constexpr index_t N{detray::traits::max_rank<M>};

  auto L = zero<M>();

  // Cholesky–Banachiewicz algorithm
  for (std::size_t i = 0u; i < N; i++) {
    for (std::size_t j = 0u; j <= i; j++) {
      scalar_t sum = 0.f;
      for (std::size_t k = 0u; k < j; k++)
        sum += elem(L, i, k) * elem(L, j, k);

      if (i == j) {
        elem(L, i, j) = static_cast<scalar_t>(
            detray::algebra::math::sqrt(elem(m, i, i) - sum));
      } else {
        elem(L, i, j) = (1.f / elem(L, j, j) * (elem(m, i, j) - sum));
      }
    }
  }

  return L;
}

}  // namespace detray::algebra::generic::math
