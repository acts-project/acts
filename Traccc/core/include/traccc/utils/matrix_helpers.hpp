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

// System include(s)
#include <cassert>

namespace traccc {

/// @brief Compute a masked inverse of a 2x2 matrix.
///
/// For 1D measurements only the [0, 0] element of the input matrix is
/// meaningful, so we compute 1 / M[0,0] and zero out the rest. For 2D
/// measurements a regular 2x2 inverse is computed.
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

}  // namespace traccc
