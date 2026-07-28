/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// traccc include
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/container.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/trigonometric_helpers.hpp"

// detray include(s).
#include <detray/tracks/tracks.hpp>

namespace traccc {

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using free_track_parameters = detray::free_track_parameters<algebra_t>;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using bound_track_parameters = detray::bound_track_parameters<algebra_t>;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using bound_parameters_vector = detray::bound_parameters_vector<algebra_t>;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using free_vector = typename free_track_parameters<algebra_t>::vector_type;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using bound_vector = typename bound_parameters_vector<algebra_t>::vector_type;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using bound_covariance =
    typename bound_track_parameters<algebra_t>::covariance_type;

template <detray::concepts::algebra algebra_t = traccc::default_algebra>
using bound_matrix = detray::bound_matrix<algebra_t>;

/// Declare all track_parameters collection types
using bound_track_parameters_collection_types =
    collection_types<bound_track_parameters<>>;

// Wrap the phi of a track parameter vector to [-pi,pi]
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE constexpr bool normalize_angles(
    bound_vector<algebra_t>& vec) {
  traccc::scalar phi;
  traccc::scalar theta;

  const traccc::scalar in_theta{getter::element(vec, e_bound_theta, 0)};
  if (math::fmod(in_theta, 2.f * constant<traccc::scalar>::pi) == 0.f ||
      in_theta == 0.f) {
    TRACCC_WARNING_HOST_DEVICE("Hit theta pole before normalization: %f",
                               in_theta);
    return false;
  }

  std::tie(phi, theta) =
      detail::wrap_phi_theta(getter::element(vec, e_bound_phi, 0), in_theta);

  if (theta <= 0.f || theta >= 2.f * constant<traccc::scalar>::pi) {
    TRACCC_WARNING_HOST_DEVICE("Hit theta pole after normalization: %f", theta);
    return false;
  }

  // Assertions of the detray bound track parameters
  assert(math::fabs(phi) <= constant<traccc::scalar>::pi);
  assert(theta <= constant<traccc::scalar>::pi);

  getter::element(vec, e_bound_phi, 0) = phi;
  getter::element(vec, e_bound_theta, 0) = theta;

  return true;
}

// Wrap the phi of bound track parameters to [-pi,pi]
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE constexpr bool normalize_angles(
    bound_parameters_vector<algebra_t>& param_vec) {
  return normalize_angles<algebra_t>(param_vec.vector());
}

/// Covariance inflation used for track fitting
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline void inflate_covariance(
    bound_track_parameters<algebra_t>& param, const traccc::scalar inf_fac) {
  auto& cov = param.covariance();
  for (unsigned int i = 0; i < e_bound_size; i++) {
    for (unsigned int j = 0; j < e_bound_size; j++) {
      if (i == j) {
        getter::element(cov, i, i) *= inf_fac;
      } else {
        getter::element(cov, i, j) = 0.f;
      }
    }
  }
}

}  // namespace traccc
