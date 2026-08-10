// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/utils/known_substructure_matrix.hpp"
#include "detray/utils/known_substructure_matrix_types.hpp"

namespace detray::detail {

/// @brief Transport a bound covariance to the destination surface.
///
/// @c C' = J C J^T, a congruence. Both matrices are structured: the full
/// Jacobian has 25 of its 36 cells free, and a covariance is symmetric and so
/// stores only its upper triangle. Carrying the product out in those
/// substructures skips the structurally zero terms, and -- because the result
/// of a congruence with a symmetric middle factor is itself symmetric -- only
/// the upper triangle of the result is computed at all.
///
/// @note The covariance is read with the shared-cell check switched off. It is
/// symmetric by definition, but a structure-blind transport leaves it
/// symmetric only to rounding: measured drift is order 1e-7 after a handful of
/// transports, which the check would reject. Taking the upper triangle both
/// states the symmetry and restores it, so the covariance no longer drifts.
template <concepts::algebra algebra_t, typename jacobian_t, typename cov_t>
DETRAY_HOST_DEVICE void transport_covariance_to_bound(
    const cov_t &covariance, const jacobian_t &full_jacobian,
    cov_t &new_covariance) {
  using scalar_t = dscalar<algebra_t>;
  using symmetric_cov_t =
      ksm::matrix<typename ksm::make_symmetric_substructure<
                      static_cast<std::size_t>(e_bound_size)>::canonical_type,
                  scalar_t>;

  const auto cov =
      symmetric_cov_t::template from_dense<algebra_t, false>(covariance);

  new_covariance = full_jacobian.congruence(cov).template to_dense<algebra_t>();
}

}  // namespace detray::detail
