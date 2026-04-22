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

// Project include(s)
#include "detray/test/common/track_generators.hpp"

#include "detray/definitions/units.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/tracks/trajectories.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/statistics.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using vector3 = test::vector3;
using point3 = test::point3;

GTEST_TEST(detray_simulation, uniform_track_generator) {
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  constexpr const scalar tol{1e-5f};
  constexpr const scalar max_pi{generator_t::configuration::k_max_pi};

  constexpr std::size_t phi_steps{50u};
  constexpr std::size_t theta_steps{50u};
  constexpr const scalar charge{-1.f};

  std::array<vector3, phi_steps * theta_steps> momenta{};

  // Loop over theta values [0,pi)
  for (std::size_t itheta{0u}; itheta < theta_steps; ++itheta) {
    const scalar theta{static_cast<scalar>(itheta) * max_pi /
                       static_cast<scalar>(theta_steps - 1u)};

    // Loop over phi values [-pi, pi)
    for (std::size_t iphi{0u}; iphi < phi_steps; ++iphi) {
      // The direction
      const scalar phi{-constant<scalar>::pi +
                       static_cast<scalar>(iphi) *
                           (constant<scalar>::pi + max_pi) /
                           static_cast<scalar>(phi_steps)};

      // initialize a track
      vector3 mom{std::cos(phi) * std::sin(theta),
                  std::sin(phi) * std::sin(theta), std::cos(theta)};
      mom = vector::normalize(mom);
      free_track_parameters<test_algebra> traj({0.f, 0.f, 0.f}, 0.f, mom, -1.f);

      momenta[itheta * phi_steps + iphi] = traj.mom(-1.f);
    }
  }

  // Now run the track generator and compare
  std::size_t n_tracks{0u};
  generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.charge(charge);
  trk_gen_cfg.phi_steps(phi_steps).theta_steps(theta_steps);

  for (const auto track : generator_t(trk_gen_cfg)) {
    vector3& expected = momenta[n_tracks];
    vector3 result = track.mom(trk_gen_cfg.charge());

    // Compare with for loop
    EXPECT_NEAR(vector::norm(expected - result), 0.f, tol)
        << "Track: \n"
        << expected[0] << "\t" << result[0] << "\n"
        << expected[1] << "\t" << result[1] << "\n"
        << expected[2] << "\t" << result[2] << std::endl;

    ++n_tracks;
  }
  ASSERT_EQ(momenta.size(), n_tracks);

  // Generate rays
  n_tracks = 0u;
  for (const auto r : generator_t(phi_steps, theta_steps)) {
    vector3& expected = momenta[n_tracks];
    vector3 result = r.dir();

    // Compare with for loop
    EXPECT_NEAR(vector::norm(expected - result), 0.f, tol)
        << "Ray: \n"
        << expected[0] << "\t" << result[0] << "\n"
        << expected[1] << "\t" << result[1] << "\n"
        << expected[2] << "\t" << result[2] << std::endl;

    ++n_tracks;
  }
  ASSERT_EQ(momenta.size(), n_tracks);

  // Generate helical trajectories
  const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                  2.f * unit<scalar>::T};
  n_tracks = 0u;
  for (const auto track : generator_t(phi_steps, theta_steps)) {
    detail::helix<test_algebra> helix_traj(track, B);
    vector3& expected = momenta[n_tracks];
    vector3 result = helix_traj.dir(0.f);

    // Compare with for loop
    EXPECT_NEAR(vector::norm(expected - result), 0.f, tol)
        << "Helix: \n"
        << expected[0] << "\t" << result[0] << "\n"
        << expected[1] << "\t" << result[1] << "\n"
        << expected[2] << "\t" << result[2] << std::endl;

    ++n_tracks;
  }
  ASSERT_EQ(momenta.size(), n_tracks);
}

GTEST_TEST(detray_simulation, uniform_track_generator_eta) {
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  constexpr const scalar tol{1e-5f};
  constexpr const scalar max_pi{generator_t::configuration::k_max_pi};
  constexpr const scalar charge{-1.f};

  constexpr std::size_t phi_steps{50u};
  constexpr std::size_t eta_steps{50u};

  std::array<vector3, phi_steps * eta_steps> momenta{};

  // Loop over eta values [-5, 5]
  for (std::size_t ieta{0u}; ieta < eta_steps; ++ieta) {
    const scalar eta{-5.f + static_cast<scalar>(ieta) * 10.f /
                                static_cast<scalar>(eta_steps - 1u)};
    const scalar theta{2.f * std::atan(std::exp(-eta))};

    // Loop over phi values [-pi, pi)
    for (std::size_t iphi{0u}; iphi < phi_steps; ++iphi) {
      // The direction
      const scalar phi{-constant<scalar>::pi +
                       static_cast<scalar>(iphi) *
                           (constant<scalar>::pi + max_pi) /
                           static_cast<scalar>(phi_steps)};

      // initialize a track
      vector3 mom{std::cos(phi) * std::sin(theta),
                  std::sin(phi) * std::sin(theta), std::cos(theta)};
      mom = vector::normalize(mom);
      free_track_parameters<test_algebra> traj({0.f, 0.f, 0.f}, 0.f, mom, -1.f);

      momenta[ieta * phi_steps + iphi] = traj.mom(charge);
    }
  }

  // Now run the track generator and compare
  std::size_t n_tracks{0u};

  auto trk_generator = generator_t{};
  trk_generator.config().phi_steps(phi_steps).eta_steps(eta_steps);
  trk_generator.config().charge(charge);

  for (const auto track : trk_generator) {
    vector3& expected = momenta[n_tracks];
    vector3 result = track.mom(trk_generator.config().charge());

    // Compare with for loop
    EXPECT_NEAR(vector::norm(expected - result), 0.f, tol)
        << "Expected\tResult: \n"
        << expected[0] << "\t" << result[0] << "\n"
        << expected[1] << "\t" << result[1] << "\n"
        << expected[2] << "\t" << result[2] << std::endl;

    ++n_tracks;
  }
  ASSERT_EQ(momenta.size(), n_tracks);
}

GTEST_TEST(detray_simulation, uniform_track_generator_with_range) {
  using generator_t =
      uniform_track_generator<free_track_parameters<test_algebra>>;

  constexpr const scalar tol{1e-5f};

  std::vector<std::array<scalar, 2>> theta_phi;

  auto trk_gen_cfg = generator_t::configuration{};
  trk_gen_cfg.phi_range(-2.f, 2.f).phi_steps(4u);
  trk_gen_cfg.theta_range(1.f, 2.f).theta_steps(2u);

  for (const auto track : generator_t{trk_gen_cfg}) {
    const auto dir = track.dir();
    theta_phi.push_back({vector::theta(dir), vector::phi(dir)});
  }

  EXPECT_EQ(theta_phi.size(), 8u);
  EXPECT_NEAR(theta_phi[0][0], 1.f, tol);
  EXPECT_NEAR(theta_phi[0][1], -2.f, tol);
  EXPECT_NEAR(theta_phi[1][0], 1.f, tol);
  EXPECT_NEAR(theta_phi[1][1], -1.f, tol);
  EXPECT_NEAR(theta_phi[2][0], 1.f, tol);
  EXPECT_NEAR(theta_phi[2][1], 0.f, tol);
  EXPECT_NEAR(theta_phi[3][0], 1.f, tol);
  EXPECT_NEAR(theta_phi[3][1], 1.f, tol);
  EXPECT_NEAR(theta_phi[4][0], 2.f, tol);
  EXPECT_NEAR(theta_phi[4][1], -2.f, tol);
  EXPECT_NEAR(theta_phi[5][0], 2.f, tol);
  EXPECT_NEAR(theta_phi[5][1], -1.f, tol);
  EXPECT_NEAR(theta_phi[6][0], 2.f, tol);
  EXPECT_NEAR(theta_phi[6][1], 0.f, tol);
  EXPECT_NEAR(theta_phi[7][0], 2.f, tol);
  EXPECT_NEAR(theta_phi[7][1], 1.f, tol);
}

/// Tests a random number based track state generator - uniform distribution
GTEST_TEST(detray_simulation, random_track_generator_uniform) {
  // Use deterministic random number generator for testing
  using uniform_gen_t =
      detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
  using trk_generator_t =
      random_track_generator<free_track_parameters<test_algebra>,
                             uniform_gen_t>;

  // Tolerance depends on sample size
  constexpr scalar tol{0.02f};

  // Track counter
  std::size_t n_tracks{0u};
  constexpr std::size_t n_gen_tracks{10000u};

  // Other params
  trk_generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.n_tracks(n_gen_tracks);
  trk_gen_cfg.seed(42u);
  trk_gen_cfg.phi_range(-0.9f * constant<scalar>::pi,
                        0.8f * constant<scalar>::pi);
  trk_gen_cfg.eta_range(-4.f, 4.f);
  trk_gen_cfg.mom_range(1.f * unit<scalar>::GeV, 2.f * unit<scalar>::GeV);
  trk_gen_cfg.origin_stddev(0.1f * unit<scalar>::mm, 0.f * unit<scalar>::mm,
                            0.2f * unit<scalar>::mm);

  // Catch the results
  std::array<scalar, n_gen_tracks> x{};
  std::array<scalar, n_gen_tracks> y{};
  std::array<scalar, n_gen_tracks> z{};
  std::array<scalar, n_gen_tracks> mom{};
  std::array<scalar, n_gen_tracks> phi{};
  std::array<scalar, n_gen_tracks> theta{};

  for (const auto track : trk_generator_t{trk_gen_cfg}) {
    const auto pos = track.pos();
    x[n_tracks] = pos[0];
    y[n_tracks] = pos[1];
    z[n_tracks] = pos[2];
    mom[n_tracks] = track.p(trk_gen_cfg.charge());
    phi[n_tracks] = vector::phi(track.dir());
    theta[n_tracks] = vector::theta(track.dir());

    ++n_tracks;
  }

  ASSERT_EQ(n_gen_tracks, n_tracks);

  // Check uniform distribution
  const auto& ori = trk_gen_cfg.origin();
  const auto& ori_stddev = trk_gen_cfg.origin_stddev();
  const auto& phi_range = trk_gen_cfg.phi_range();
  const auto& theta_range = trk_gen_cfg.theta_range();
  ASSERT_NEAR(theta_range[0], 0.0366f, tol);
  ASSERT_NEAR(theta_range[1], 3.105f, tol);
  const auto& mom_range = trk_gen_cfg.mom_range();

  // Mean
  EXPECT_NEAR(statistics::mean(x), ori[0], tol);
  EXPECT_NEAR(statistics::mean(y), ori[1], tol);
  EXPECT_NEAR(statistics::mean(z), ori[2], tol);
  EXPECT_NEAR(statistics::mean(mom), 0.5f * (mom_range[0] + mom_range[1]), tol);
  EXPECT_NEAR(statistics::mean(phi), 0.5f * (phi_range[0] + phi_range[1]), tol);
  EXPECT_NEAR(statistics::mean(theta), 0.5f * (theta_range[0] + theta_range[1]),
              tol);

  // variance
  EXPECT_NEAR(statistics::variance(x), ori_stddev[0] * ori_stddev[0], tol);
  EXPECT_NEAR(statistics::variance(y), ori_stddev[1] * ori_stddev[1], tol);
  EXPECT_NEAR(statistics::variance(z), ori_stddev[2] * ori_stddev[2], tol);
  EXPECT_NEAR(statistics::variance(mom),
              1.0f / 12.0f * std::pow(mom_range[1] - mom_range[0], 2.f), tol);
  EXPECT_NEAR(statistics::variance(phi),
              1.0f / 12.0f * std::pow(phi_range[1] - phi_range[0], 2.f), tol);
  EXPECT_NEAR(statistics::variance(theta),
              1.0f / 12.0f * std::pow(theta_range[1] - theta_range[0], 2.f),
              tol);
}

/// Tests a random number based track state generator - normal distribution
GTEST_TEST(detray_simulation, random_track_generator_normal) {
  // Use deterministic random number generator for testing
  using normal_gen_t =
      detail::random_numbers<scalar, std::normal_distribution<scalar>>;
  using trk_generator_t =
      random_track_generator<free_track_parameters<test_algebra>, normal_gen_t>;

  // Tolerance depends on sample size
  constexpr scalar tol{0.02f};

  // Track counter
  std::size_t n_tracks{0u};
  constexpr std::size_t n_gen_tracks{10000u};

  // Other params
  trk_generator_t::configuration trk_gen_cfg{};
  trk_gen_cfg.n_tracks(n_gen_tracks);
  trk_gen_cfg.seed(42u);
  trk_gen_cfg.phi_range(-0.9f * constant<scalar>::pi,
                        0.8f * constant<scalar>::pi);
  trk_gen_cfg.eta_range(-4.f, 4.f);
  trk_gen_cfg.mom_range(1.f * unit<scalar>::GeV, 2.f * unit<scalar>::GeV);
  trk_gen_cfg.origin(0.f, 0.f, 0.f);
  trk_gen_cfg.origin_stddev(0.1f * unit<scalar>::mm, 0.5f * unit<scalar>::mm,
                            0.3f * unit<scalar>::mm);

  // Catch the results
  std::array<scalar, n_gen_tracks> x{};
  std::array<scalar, n_gen_tracks> y{};
  std::array<scalar, n_gen_tracks> z{};
  std::array<scalar, n_gen_tracks> mom{};
  std::array<scalar, n_gen_tracks> phi{};
  std::array<scalar, n_gen_tracks> theta{};

  for (const auto track : trk_generator_t{trk_gen_cfg}) {
    const auto pos = track.pos();
    x[n_tracks] = pos[0];
    y[n_tracks] = pos[1];
    z[n_tracks] = pos[2];
    mom[n_tracks] = track.p(trk_gen_cfg.charge());
    phi[n_tracks] = vector::phi(track.dir());
    theta[n_tracks] = vector::theta(track.dir());

    ++n_tracks;
  }

  ASSERT_EQ(n_gen_tracks, n_tracks);

  // check gaussian distribution - values are clamped to phi/theta range
  const auto& ori = trk_gen_cfg.origin();
  const auto& ori_stddev = trk_gen_cfg.origin_stddev();
  const auto& phi_range = trk_gen_cfg.phi_range();
  const auto& theta_range = trk_gen_cfg.theta_range();
  ASSERT_NEAR(theta_range[0], 0.0366f, tol);
  ASSERT_NEAR(theta_range[1], 3.105f, tol);
  const auto& mom_range = trk_gen_cfg.mom_range();

  // Mean
  EXPECT_NEAR(statistics::mean(x), ori[0], tol);
  EXPECT_NEAR(statistics::mean(y), ori[1], tol);
  EXPECT_NEAR(statistics::mean(z), ori[2], tol);
  EXPECT_NEAR(statistics::mean(mom),
              mom_range[0] + 0.5f * (mom_range[1] - mom_range[0]), tol);
  EXPECT_NEAR(statistics::mean(phi),
              phi_range[0] + 0.5f * (phi_range[1] - phi_range[0]), tol);
  EXPECT_NEAR(statistics::mean(theta),
              theta_range[0] + 0.5f * (theta_range[1] - theta_range[0]), tol);

  // variance
  EXPECT_NEAR(statistics::variance(x), ori_stddev[0] * ori_stddev[0], tol);
  EXPECT_NEAR(statistics::variance(y), ori_stddev[1] * ori_stddev[1], tol);
  EXPECT_NEAR(statistics::variance(z), ori_stddev[2] * ori_stddev[2], tol);
  EXPECT_NEAR(statistics::variance(mom),
              std::pow(0.5f / 3.0f * (mom_range[1] - mom_range[0]), 2.f), tol);
  EXPECT_NEAR(statistics::variance(phi),
              std::pow(0.5f / 3.0f * (phi_range[1] - phi_range[0]), 2.f), tol);
  EXPECT_NEAR(statistics::variance(theta),
              std::pow(0.5f / 3.0f * (theta_range[1] - theta_range[0]), 2.f),
              tol);
}
