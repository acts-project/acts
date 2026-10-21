// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/common/event_generator/random_numbers.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <ostream>
#include <random>

namespace detray {

/// Configuration for the random track generator
template <concepts::scalar scalar_t>
struct random_track_generator_config {
  using seed_t = std::uint64_t;

  /// Gaussian vertex smearing
  bool m_do_vtx_smearing = false;

  /// Monte-Carlo seed
  seed_t m_seed{detail::random_numbers<scalar_t>::default_seed()};

  /// How many tracks will be generated
  std::size_t m_n_tracks{10u};

  /// Range for phi [-pi, pi) and theta [0, pi)
  darray<scalar_t, 2> m_phi_range{-constant<scalar_t>::pi,
                                  constant<scalar_t>::pi};
  darray<scalar_t, 2> m_theta_range{0.f, constant<scalar_t>::pi};

  /// Momentum range
  darray<scalar_t, 2> m_mom_range{1.f * unit<scalar_t>::GeV,
                                  1.f * unit<scalar_t>::GeV};
  /// Whether to interpret the momentum @c m_mom_range as p_T
  bool m_is_pT{false};

  /// Track origin
  darray<scalar_t, 3> m_origin{0.f, 0.f, 0.f};
  darray<scalar_t, 3> m_origin_stddev{0.f, 0.f, 0.f};

  /// Randomly flip the charge sign?
  bool m_randomize_charge{false};

  /// Time parameter and charge of the track
  scalar_t m_time{0.f * unit<scalar_t>::us};
  scalar_t m_charge{-1.f * unit<scalar_t>::e};

  /// Setters
  /// @{
  DETRAY_HOST_DEVICE random_track_generator_config& seed(const seed_t s) {
    m_seed = s;
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& do_vertex_smearing(bool b) {
    m_do_vtx_smearing = b;
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& n_tracks(std::size_t n) {
    assert(n > 0);
    m_n_tracks = n;
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& phi_range(
      const scalar_t low, const scalar_t high) {
    auto min_phi{
        std::clamp(low, -constant<scalar_t>::pi, constant<scalar_t>::pi)};
    auto max_phi{
        std::clamp(high, -constant<scalar_t>::pi, constant<scalar_t>::pi)};

    assert(min_phi <= max_phi);

    m_phi_range = {min_phi, max_phi};
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE random_track_generator_config& phi_range(
      darray<o_scalar_t, 2> r) {
    phi_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& theta_range(scalar_t low,
                                                                scalar_t high) {
    auto min_theta{std::clamp(low, scalar_t{0.f}, constant<scalar_t>::pi)};
    auto max_theta{std::clamp(high, scalar_t{0.f}, constant<scalar_t>::pi)};

    assert(min_theta <= max_theta);

    m_theta_range = {min_theta, max_theta};
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE random_track_generator_config& theta_range(
      darray<o_scalar_t, 2> r) {
    theta_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& eta_range(scalar_t low,
                                                              scalar_t high) {
    // This value is more or less random
    constexpr auto num_max{0.001f * std::numeric_limits<scalar_t>::max()};
    auto min_eta{low > -num_max ? low : -num_max};
    auto max_eta{high < num_max ? high : num_max};

    assert(min_eta <= max_eta);

    auto get_theta = [](const scalar_t eta) {
      return 2.f * math::atan(math::exp(-eta));
    };

    theta_range(get_theta(max_eta), get_theta(min_eta));
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE random_track_generator_config& eta_range(
      darray<o_scalar_t, 2> r) {
    eta_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& mom_range(scalar_t low,
                                                              scalar_t high) {
    m_is_pT = false;
    assert(low >= 0.f);
    assert(low <= high);
    m_mom_range = {low, high};
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE random_track_generator_config& mom_range(
      darray<o_scalar_t, 2> r) {
    mom_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& pT_range(scalar_t low,
                                                             scalar_t high) {
    m_is_pT = true;
    assert(low >= 0.f);
    assert(low <= high);
    m_mom_range = {low, high};
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE random_track_generator_config& pT_range(
      darray<o_scalar_t, 2> r) {
    pT_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& p_tot(scalar_t p) {
    mom_range(p, p);
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& p_T(scalar_t p) {
    pT_range(p, p);
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& origin(const scalar_t x,
                                                           const scalar_t y,
                                                           const scalar_t z) {
    m_origin = {x, y, z};
    return *this;
  }
  template <concepts::point3D point3_t>
  DETRAY_HOST_DEVICE random_track_generator_config& origin(
      const point3_t& ori) {
    return origin(ori[0], ori[1], ori[2]);
  }
  DETRAY_HOST_DEVICE random_track_generator_config& origin_stddev(
      const scalar_t x, const scalar_t y, const scalar_t z) {
    m_do_vtx_smearing = true;
    m_origin_stddev = {x, y, z};
    return *this;
  }
  template <concepts::point3D point3_t>
  DETRAY_HOST_DEVICE random_track_generator_config& origin_stddev(
      point3_t stddev) {
    return origin_stddev(stddev[0], stddev[1], stddev[2]);
  }
  DETRAY_HOST_DEVICE
  random_track_generator_config& randomize_charge(bool rc) {
    m_randomize_charge = rc;
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& time(scalar_t t) {
    assert(t >= 0.f);
    m_time = t;
    return *this;
  }
  DETRAY_HOST_DEVICE random_track_generator_config& charge(scalar_t q) {
    m_charge = q;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  DETRAY_HOST_DEVICE constexpr seed_t seed() const { return m_seed; }
  DETRAY_HOST_DEVICE constexpr bool do_vertex_smearing() const {
    return m_do_vtx_smearing;
  }
  DETRAY_HOST_DEVICE constexpr std::size_t n_tracks() const {
    return m_n_tracks;
  }
  DETRAY_HOST_DEVICE constexpr const darray<scalar_t, 2>& phi_range() const {
    return m_phi_range;
  }
  DETRAY_HOST_DEVICE constexpr const darray<scalar_t, 2>& theta_range() const {
    return m_theta_range;
  }
  DETRAY_HOST_DEVICE constexpr const darray<scalar_t, 2>& mom_range() const {
    return m_mom_range;
  }
  DETRAY_HOST_DEVICE constexpr const auto& origin() const { return m_origin; }
  DETRAY_HOST_DEVICE constexpr const auto& origin_stddev() const {
    return m_origin_stddev;
  }
  DETRAY_HOST_DEVICE constexpr bool is_pT() const { return m_is_pT; }
  DETRAY_HOST_DEVICE constexpr bool randomize_charge() const {
    return m_randomize_charge;
  }
  DETRAY_HOST_DEVICE constexpr scalar_t time() const { return m_time; }
  DETRAY_HOST_DEVICE constexpr scalar_t charge() const { return m_charge; }
  /// @}

  /// Print the random track generator configuration
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out,
                                  const random_track_generator_config& cfg) {
    const auto& ori = cfg.origin();
    const auto& mom_range = cfg.mom_range();
    const auto& phi_range = cfg.phi_range();
    const auto& theta_range = cfg.theta_range();

    // General
    out << "\nRandom track generator\n"
        << "----------------------------\n"
        << "  Random seed           : " << cfg.seed() << "\n"
        << "  No. tracks            : " << cfg.n_tracks() << "\n"
        << "  Charge                : "
        << cfg.charge() / detray::unit<scalar_t>::e << " [e]\n"
        << "  Rand. charge          : " << std::boolalpha
        << cfg.randomize_charge() << std::noboolalpha << "\n";

    // Momentum and direction
    if (cfg.is_pT()) {
      out << "  Transverse mom.       : [";
    } else {
      out << "  Momentum              : [";
    }
    out << mom_range[0] / detray::unit<scalar_t>::GeV << ", "
        << mom_range[1] / detray::unit<scalar_t>::GeV << ") [GeV]\n"
        << "  Phi range             : ["
        << phi_range[0] / detray::unit<scalar_t>::rad << ", "
        << phi_range[1] / detray::unit<scalar_t>::rad << ") [rad]\n"
        << "  Theta range           : ["
        << theta_range[0] / detray::unit<scalar_t>::rad << ", "
        << theta_range[1] / detray::unit<scalar_t>::rad << ") [rad]\n"
        << "  Origin                : [" << ori[0] / detray::unit<scalar_t>::mm
        << ", " << ori[1] / detray::unit<scalar_t>::mm << ", "
        << ori[2] / detray::unit<scalar_t>::mm << "] [mm]\n"
        << "  Do vertex smearing    : " << std::boolalpha
        << cfg.do_vertex_smearing() << "\n"
        << std::noboolalpha;

    if (cfg.do_vertex_smearing()) {
      const auto& ori_stddev = cfg.origin_stddev();
      out << "  Origin stddev         : ["
          << ori_stddev[0] / detray::unit<scalar_t>::mm << ", "
          << ori_stddev[1] / detray::unit<scalar_t>::mm << ", "
          << ori_stddev[2] / detray::unit<scalar_t>::mm << "] [mm]\n";
    }

    return out;
  }
};

}  // namespace detray
