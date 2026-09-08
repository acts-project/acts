/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/utils/concepts.hpp"

// Detray include(s).
#include <detray/algebra/common/known_substructure_matrix.hpp>

// System include(s)
#include <cassert>
#include <cstddef>
#include <utility>

namespace traccc {

/// @brief Compute a masked inverse of a dense 2x2 matrix.
///
/// For 1D measurements only the [0, 0] element of the input matrix is
/// meaningful, so we compute 1 / M[0,0] and zero out the rest. For 2D
/// measurements a regular 2x2 inverse is computed.
///
/// @note This is the dense spelling, kept for the code paths that do not run
/// on structured matrices yet. Prefer the overload below, which needs no
/// conversion on either side.
///
/// @param M   the 2x2 matrix to invert
/// @param dim the meaningful dimension of @c M (1 or 2)
/// @returns the masked inverse of @c M
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline detray::dmatrix<algebra_t, 2, 2> masked_inverse(
    const detray::dmatrix<algebra_t, 2, 2>& M, const unsigned int dim) {
  assert(dim == 1u || dim == 2u);

  detray::dmatrix<algebra_t, 2, 2> M_inv;
  if (dim == 1u) {
    assert(getter::element(M, 0u, 0u) != 0.f);
    getter::element(M_inv, 0u, 0u) = 1.f / getter::element(M, 0u, 0u);
    getter::element(M_inv, 0u, 1u) = 0.f;
    getter::element(M_inv, 1u, 0u) = 0.f;
    getter::element(M_inv, 1u, 1u) = 0.f;
  } else {
    M_inv = matrix::inverse(M);
  }
  return M_inv;
}

/// @brief The substructures of the matrices a Kalman update runs on
///
/// These are the three facts that make the products of an update cheap. The
/// library skips the multiplications by a structural zero, so naming the
/// structure here is what removes the work.
/// @{
/// A track covariance, which is symmetric.
template <std::size_t N>
using track_covariance_substructure =
    typename detray::ksm::make_symmetric_substructure<N>::canonical_type;

/// An observation model. A measurement subspace only ever selects
/// @c e_bound_loc0 and @c e_bound_loc1, so every further column is zero.
template <std::size_t D, std::size_t N>
using observation_model_substructure =
    typename detray::ksm::make_left_columns_substructure<D, N,
                                                         2u>::canonical_type;

/// A measurement covariance, which is diagonal:
/// @c edm::get_measurement_covariance writes its off-diagonal elements as
/// exact zeros.
template <std::size_t D>
using measurement_covariance_substructure =
    typename detray::ksm::make_diagonal_substructure<D>::canonical_type;

/// A residual-space covariance and its inverse, which is symmetric but not
/// diagonal: R = V - H*C*H^T is full even though V is diagonal, so the
/// measurement covariance shape would drop the off-diagonal element.
template <std::size_t D>
using residual_covariance_substructure =
    typename detray::ksm::make_symmetric_substructure<D>::canonical_type;
/// @}

/// @brief Read a dense square matrix into a matrix of symmetric substructure.
///
/// A symmetric substructure keeps a single scalar for the pair (i, j) and
/// (j, i), so only the upper triangle of @p m is read and the lower one is
/// discarded.
///
/// @c ksm::matrix::from_dense cannot serve here, because it asserts that the
/// two mirrored elements hold the same value. A transported track covariance
/// is symmetric in exact arithmetic only: the generated
/// @c transport_covariance_to_bound_impl builds every element of J*C*J^T as
/// its own expression, so two mirrored elements differ in the last bits.
///
/// @tparam ksm_matrix_t the destination type, whose substructure must be
///                      square and symmetric
/// @param M the dense matrix to read
///
/// @returns @c M with its lower triangle discarded
template <typename ksm_matrix_t, typename dense_t>
TRACCC_HOST_DEVICE inline ksm_matrix_t symmetric_from_dense(const dense_t& M) {
  using substructure_t = typename ksm_matrix_t::canonical_substructure_type;
  using scalar_t = typename ksm_matrix_t::scalar_type;

  static_assert(detray::ksm::concepts::is_symmetric<substructure_t>,
                "the destination substructure is not symmetric");

  constexpr std::size_t N = substructure_t::rows;
  static_assert(static_cast<std::size_t>(detray::traits::rows<dense_t>) == N);
  static_assert(static_cast<std::size_t>(detray::traits::columns<dense_t>) ==
                N);

  ksm_matrix_t rv{};

  const auto read = [&rv, &M]<std::size_t I, std::size_t J>() {
    if constexpr (I <= J) {
      rv.template at<I, J>() = static_cast<scalar_t>(getter::element<I, J>(M));
    }
  };

  [&read]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (read.template operator()<Ks / N, Ks % N>(), ...);
  }(std::make_index_sequence<N * N>{});

  return rv;
}

/// @brief Compute a masked inverse of a symmetric 2x2 matrix.
///
/// Only the leading @c DIM by @c DIM block of @p m is meaningful, and the rest
/// of the result is zero. For a 1D measurement that leaves a single
/// reciprocal; the second row and column of the input hold a padding variance
/// that must not reach the result.
///
/// The dimension is a template parameter so that each case is its own body:
/// the 1D case is one division, and writing it as a runtime mask would make
/// it pay for the 2x2 inverse it does not need.
///
/// @tparam DIM the meaningful dimension of @p m, 1 or 2
/// @param m the matrix to invert
/// @returns the masked inverse of @p m
template <std::size_t DIM, concepts::symmetric_matrix<2> matrix_t>
TRACCC_HOST_DEVICE constexpr matrix_t masked_inverse(const matrix_t& m) {
  static_assert(DIM == 1u || DIM == 2u, "a measurement is 1D or 2D");

  // @c concepts::symmetric_matrix is only a readability hint for a
  // @c ksm::matrix, so the substructure is what has to be asked. A dense 2x2
  // would pass the concept and then leave [1, 0] unwritten below.
  static_assert(detray::ksm::concepts::is_symmetric<
                    typename matrix_t::canonical_substructure_type>,
                "the matrix to invert is not of a symmetric substructure");

  if constexpr (DIM == 2u) {
    // The whole matrix is meaningful, so this is an ordinary inverse. The
    // result is exactly symmetric, because the two mirrored cells are one
    // stored variable.
    return m.inverse();
  } else {
    // Doubled parentheses: the comma in the template argument list would
    // otherwise be read as a macro argument separator.
    assert((m.template at<0, 0>() != 0.f));

    matrix_t rv{};
    rv.template at<0, 0>() = 1.f / m.template at<0, 0>();
    rv.template at<0, 1>() = 0.f;
    rv.template at<1, 1>() = 0.f;
    return rv;
  }
}

/// @brief Compute a masked inverse of a symmetric 2x2 matrix.
///
/// The measurement dimension is known only at run time, so this selects
/// between the two bodies above.
///
/// @param m the matrix to invert
/// @param dim the meaningful dimension of @p m (1 or 2)
/// @returns the masked inverse of @p m
template <concepts::symmetric_matrix<2> matrix_t>
TRACCC_HOST_DEVICE constexpr matrix_t masked_inverse(const matrix_t& m,
                                                     const unsigned int dim) {
  assert(dim == 1u || dim == 2u);

  return (dim == 1u) ? masked_inverse<1u>(m) : masked_inverse<2u>(m);
}

}  // namespace traccc
