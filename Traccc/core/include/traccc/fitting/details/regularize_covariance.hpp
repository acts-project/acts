// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/track_parameters.hpp"
#include "traccc/utils/logging.hpp"

namespace traccc::details {

/// Check the covariance matrix and try to make it positive semi-definite
///
/// @param[out] cov  covariance matrix
/// @param[in] min_var variance threshold below which to flag an error
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE constexpr bool regularize_covariance(
    traccc::bound_matrix<algebra_t>& cov,
    const detray::dscalar<algebra_t> min_var) {
  if (getter::element(cov, 0, 0) < min_var ||
      getter::element(cov, 1, 1) < min_var ||
      getter::element(cov, 2, 2) < min_var ||
      getter::element(cov, 3, 3) < min_var ||
      getter::element(cov, 4, 4) < min_var ||
      getter::element(cov, 5, 5) < min_var) {
    TRACCC_ERROR_HOST_DEVICE("Negative variance");
    return false;
  } else if (getter::element(cov, 0, 0) < 0.f ||
             getter::element(cov, 1, 1) < 0.f ||
             getter::element(cov, 2, 2) < 0.f ||
             getter::element(cov, 3, 3) < 0.f ||
             getter::element(cov, 4, 4) < 0.f ||
             getter::element(cov, 5, 5) < 0.f) {
    TRACCC_WARNING_HOST_DEVICE("Negative variance: Regularize...");
  }

  getter::element(cov, 0, 0) = math::fabs(getter::element(cov, 0, 0));
  getter::element(cov, 1, 1) = math::fabs(getter::element(cov, 1, 1));
  getter::element(cov, 2, 2) = math::fabs(getter::element(cov, 2, 2));
  getter::element(cov, 3, 3) = math::fabs(getter::element(cov, 3, 3));
  getter::element(cov, 4, 4) = math::fabs(getter::element(cov, 4, 4));
  getter::element(cov, 5, 5) = math::fabs(getter::element(cov, 5, 5));

  return true;
}

}  // namespace traccc::details
