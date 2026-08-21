/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/type_traits.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/utils/concepts.hpp"
#include "traccc/utils/logging.hpp"

namespace traccc::details {

/// Check the covariance matirx and try to make it positive semi-definite
///
/// @param[out] cov  covariance matrix
/// @param[in] min_var variance threshold below which to flag an error
template <concepts::symmetric_matrix<e_bound_size> matrix_t>
TRACCC_HOST_DEVICE constexpr bool regularize_covariance(
    matrix_t& cov, const detray::traits::scalar_t<matrix_t> min_var) {
  if (getter::element<0, 0>(cov) < min_var ||
      getter::element<1, 1>(cov) < min_var ||
      getter::element<2, 2>(cov) < min_var ||
      getter::element<3, 3>(cov) < min_var ||
      getter::element<4, 4>(cov) < min_var ||
      getter::element<5, 5>(cov) < min_var) {
    TRACCC_ERROR_HOST_DEVICE("Negative variance");
    return false;
  } else if (getter::element<0, 0>(cov) < 0.f ||
             getter::element<1, 1>(cov) < 0.f ||
             getter::element<2, 2>(cov) < 0.f ||
             getter::element<3, 3>(cov) < 0.f ||
             getter::element<4, 4>(cov) < 0.f ||
             getter::element<5, 5>(cov) < 0.f) {
    TRACCC_WARNING_HOST_DEVICE("Negative variance: Regularize...");
  }

  getter::element<0, 0>(cov) = math::fabs(getter::element<0, 0>(cov));
  getter::element<1, 1>(cov) = math::fabs(getter::element<1, 1>(cov));
  getter::element<2, 2>(cov) = math::fabs(getter::element<2, 2>(cov));
  getter::element<3, 3>(cov) = math::fabs(getter::element<3, 3>(cov));
  getter::element<4, 4>(cov) = math::fabs(getter::element<4, 4>(cov));
  getter::element<5, 5>(cov) = math::fabs(getter::element<5, 5>(cov));

  return true;
}

}  // namespace traccc::details
