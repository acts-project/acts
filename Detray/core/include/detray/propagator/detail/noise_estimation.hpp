// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/logging.hpp"

namespace detray::detail {

/// Estimate the mask tolerance that is needed for the navigation to include
/// the next surface from the current covariance
template <concepts::algebra algebra_t, typename propagator_state_t>
DETRAY_HOST_DEVICE constexpr void estimate_external_mask_tolerance(
    const bound_track_parameters<algebra_t>& bound_params,
    propagator_state_t& propagation, const dscalar<algebra_t> n_stddev,
    const dscalar<algebra_t> accumulated_error = 0.f) {
  DETRAY_VERBOSE_HOST_DEVICE("Estimate noise due to scattering...");

  using scalar_t = dscalar<algebra_t>;

  auto& navigation = propagation.navigation();
  auto& stepping = propagation.stepping();

  DETRAY_VERBOSE_HOST_DEVICE(
      "-> Is material surface: %s",
      (navigation.encountered_sf_material() ? "yes" : "no"));

  // Set the noise to be expected by the navigator after all of the actors
  // are done and the covariance is up to date

  // Positional error on the current surface as an estimate of the error
  // on the next surface
  const auto& cov = bound_params.covariance();

  const scalar_t var_loc0{getter::element(cov, e_bound_loc0, e_bound_loc0)};
  const scalar_t var_loc1{getter::element(cov, e_bound_loc1, e_bound_loc1)};

  DETRAY_DEBUG_HOST_DEVICE("-> delta loc0: %f",
                           n_stddev * math::sqrt(var_loc0));
  DETRAY_DEBUG_HOST_DEVICE("-> delta loc1: %f",
                           n_stddev * math::sqrt(var_loc1));

  // Rough estimation of the track displacement at the next surface
  const scalar_t delta_phi{
      n_stddev * math::sqrt(getter::element(cov, e_bound_phi, e_bound_phi))};
  const scalar_t delta_theta{
      n_stddev *
      math::sqrt(getter::element(cov, e_bound_theta, e_bound_theta))};

  DETRAY_DEBUG_HOST_DEVICE("-> delta phi: %f", delta_phi);
  DETRAY_DEBUG_HOST_DEVICE("-> delta theta: %f", delta_theta);

  // Calculate the difference in cartesian coordinates, as the conversion
  // uses less trigonometric functions than calculating the distance in
  // spherical coordinates
  const scalar_t phi_err{bound_params.phi() + delta_phi};
  const scalar_t theta_err{bound_params.theta() + delta_theta};
  const scalar_t sin_theta_err{math::sin(theta_err)};

  dvector3D<algebra_t> displ{math::cos(phi_err) * sin_theta_err,
                             math::sin(phi_err) * sin_theta_err,
                             math::cos(theta_err)};

  // Guess the portal envelope distance if there is no next target
  constexpr auto max_tol{5.f * unit<scalar_t>::mm};
  const scalar_t path{
      navigation.cache_exhausted()
          ? max_tol
          : math::fabs(std::as_const(navigation).target().path())};

  displ = path * (displ - stepping().dir());

  DETRAY_DEBUG_HOST_DEVICE("-> estimated displacement: %f",
                           math::sqrt(vector::dot(displ, displ)));

  // Parametrized noise component that scales with the path length
  // Accounts for material/interaction mismodelling
  const scalar_t q{stepping.particle_hypothesis().charge()};
  const scalar_t accumulated_noise{1.f * unit<scalar_t>::GeV / stepping().p(q) *
                                   accumulated_error * stepping.path_length()};

  DETRAY_DEBUG_HOST_DEVICE("-> accumulated noise: %f", accumulated_noise);

  const scalar_t ext_tol{math::sqrt(
      n_stddev * n_stddev * (var_loc0 + var_loc1) + vector::dot(displ, displ) +
      accumulated_noise * accumulated_noise)};

  // Clip to 5mm if the covariances are very large
  navigation.set_external_tol(ext_tol > max_tol ? max_tol : ext_tol);

  assert(std::isfinite(navigation.external_tol()));

  DETRAY_DEBUG_HOST_DEVICE("=> scattering noise: %f",
                           navigation.external_tol());
}

}  // namespace detray::detail
