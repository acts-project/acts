// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/common/event_generator/random_numbers.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <ostream>

namespace detray {

/// Configuration for the uniform track generator
template <concepts::scalar scalar_t>
struct uniform_track_generator_config {
  using seed_t = std::uint64_t;

  /// Monte-Carlo seed
  seed_t m_seed{detail::random_numbers<scalar_t>::default_seed()};

  /// Ensure same angle space as random track generator
  static constexpr scalar_t k_max_pi{constant<scalar_t>::pi -
                                     std::numeric_limits<scalar_t>::epsilon()};

  /// Range for phi [-pi, pi) and theta [0, pi)
  darray<scalar_t, 2> m_phi_range{-constant<scalar_t>::pi, k_max_pi};
  darray<scalar_t, 2> m_theta_range{0.f, k_max_pi};
  darray<scalar_t, 2> m_eta_range{-5.f, 5.f};

  /// Angular step size
  std::size_t m_phi_steps{50u};
  std::size_t m_theta_steps{50u};

  /// Do uniform eta steps instead of uniform theta steps
  /// (use same number of steps and range)
  bool m_uniform_eta{false};

  /// Track origin
  darray<scalar_t, 3> m_origin{0.f, 0.f, 0.f};

  /// Magnitude of momentum: Default is one to keep directions normalized
  /// if no momentum information is needed (e.g. for a ray)
  scalar_t m_p_mag{1.f * unit<scalar_t>::GeV};
  /// Whether to interpret the momentum @c m_p_mag as p_T
  bool m_is_pT{false};

  /// Randomly flip the charge sign?
  bool m_randomize_charge{false};

  /// Time parameter and charge of the track
  scalar_t m_time{0.f * unit<scalar_t>::us};
  scalar_t m_charge{-1.f * unit<scalar_t>::e};

  /// Setters
  /// @{
  DETRAY_HOST_DEVICE uniform_track_generator_config& seed(const seed_t s) {
    m_seed = s;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& phi_range(scalar_t low,
                                                               scalar_t high) {
    auto min_phi{std::clamp(low, -constant<scalar_t>::pi, k_max_pi)};
    auto max_phi{std::clamp(high, -constant<scalar_t>::pi, k_max_pi)};

    assert(min_phi <= max_phi);

    m_phi_range = {min_phi, max_phi};
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE uniform_track_generator_config& phi_range(
      darray<o_scalar_t, 2> r) {
    phi_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& theta_range(
      scalar_t low, scalar_t high) {
    auto min_theta{std::clamp(low, scalar_t{0.f}, k_max_pi)};
    auto max_theta{std::clamp(high, scalar_t{0.f}, k_max_pi)};

    assert(min_theta <= max_theta);

    m_theta_range = {min_theta, max_theta};
    m_uniform_eta = false;
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE uniform_track_generator_config& theta_range(
      darray<o_scalar_t, 2> r) {
    theta_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& eta_range(scalar_t low,
                                                               scalar_t high) {
    // This value is more or less random
    constexpr auto num_max{0.001f * std::numeric_limits<scalar_t>::max()};
    auto min_eta{low > -num_max ? low : -num_max};
    auto max_eta{high < num_max ? high : num_max};

    assert(min_eta <= max_eta);

    m_eta_range = {min_eta, max_eta};
    m_uniform_eta = true;
    return *this;
  }
  template <typename o_scalar_t>
  DETRAY_HOST_DEVICE uniform_track_generator_config& eta_range(
      darray<o_scalar_t, 2> r) {
    eta_range(static_cast<o_scalar_t>(r[0]), static_cast<o_scalar_t>(r[1]));
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& phi_steps(std::size_t n) {
    assert(n > 0);
    m_phi_steps = n;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& theta_steps(
      std::size_t n) {
    assert(n > 0);
    m_theta_steps = n;
    m_uniform_eta = false;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& eta_steps(std::size_t n) {
    assert(n > 0);
    m_theta_steps = n;
    m_uniform_eta = true;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& uniform_eta(bool b) {
    m_uniform_eta = b;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& origin(const scalar_t x,
                                                            const scalar_t y,
                                                            const scalar_t z) {
    m_origin = {x, y, z};
    return *this;
  }
  template <concepts::point3D point3_t>
  DETRAY_HOST_DEVICE uniform_track_generator_config& origin(
      const point3_t& ori) {
    return origin(ori[0], ori[1], ori[2]);
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& p_tot(scalar_t p) {
    assert(p > 0.f);
    m_is_pT = false;
    m_p_mag = p;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& p_T(scalar_t p) {
    assert(p > 0.f);
    m_is_pT = true;
    m_p_mag = p;
    return *this;
  }
  DETRAY_HOST_DEVICE
  uniform_track_generator_config& randomize_charge(bool rc) {
    m_randomize_charge = rc;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& time(scalar_t t) {
    m_time = t;
    return *this;
  }
  DETRAY_HOST_DEVICE uniform_track_generator_config& charge(scalar_t q) {
    m_charge = q;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  DETRAY_HOST_DEVICE constexpr seed_t seed() const { return m_seed; }
  DETRAY_HOST_DEVICE constexpr std::size_t n_tracks() const {
    return phi_steps() * theta_steps();
  }
  DETRAY_HOST_DEVICE constexpr darray<scalar_t, 2> phi_range() const {
    return m_phi_range;
  }
  DETRAY_HOST_DEVICE constexpr darray<scalar_t, 2> theta_range() const {
    return m_theta_range;
  }
  DETRAY_HOST_DEVICE constexpr darray<scalar_t, 2> eta_range() const {
    return m_eta_range;
  }
  DETRAY_HOST_DEVICE constexpr darray<scalar_t, 2> mom_range() const {
    return {m_p_mag, m_p_mag};
  }
  DETRAY_HOST_DEVICE constexpr std::size_t phi_steps() const {
    return m_phi_steps;
  }
  DETRAY_HOST_DEVICE constexpr std::size_t theta_steps() const {
    return m_theta_steps;
  }
  DETRAY_HOST_DEVICE constexpr std::size_t eta_steps() const {
    return m_theta_steps;
  }
  DETRAY_HOST_DEVICE constexpr bool uniform_eta() const {
    return m_uniform_eta;
  }
  DETRAY_HOST_DEVICE constexpr const auto& origin() const { return m_origin; }
  DETRAY_HOST_DEVICE constexpr bool is_pT() const { return m_is_pT; }
  DETRAY_HOST_DEVICE constexpr bool randomize_charge() const {
    return m_randomize_charge;
  }
  DETRAY_HOST_DEVICE constexpr scalar_t time() const { return m_time; }
  DETRAY_HOST_DEVICE constexpr scalar_t charge() const { return m_charge; }
  /// @}

  /// Print the uniform track generator configuration
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out,
                                  const uniform_track_generator_config& cfg) {
    const auto& ori = cfg.origin();
    const auto& phi_range = cfg.phi_range();

    // General
    out << "\nUniform track generator\n"
        << "----------------------------\n"
        << "  Random seed           : " << cfg.seed() << "\n"
        << "  No. tracks            : " << cfg.n_tracks() << "\n"
        << "    -> phi steps        : " << cfg.phi_steps() << "\n"
        << "    -> theta/eta steps  : " << cfg.theta_steps() << "\n"
        << "  Charge                : "
        << cfg.charge() / detray::unit<scalar_t>::e << " [e]\n"
        << "  Rand. charge          : " << std::boolalpha
        << cfg.randomize_charge() << std::noboolalpha << "\n";

    // Momentum
    if (cfg.is_pT()) {
      out << "  Transverse mom.       : "
          << cfg.m_p_mag / detray::unit<scalar_t>::GeV << " [GeV]\n";
    } else {
      out << "  Momentum              : "
          << cfg.m_p_mag / detray::unit<scalar_t>::GeV << " [GeV]\n";
    }

    // Direction
    out << "  Phi range             : ["
        << phi_range[0] / detray::unit<scalar_t>::rad << ", "
        << phi_range[1] / detray::unit<scalar_t>::rad << ") [rad]\n";
    if (cfg.uniform_eta()) {
      const auto& eta_range = cfg.eta_range();
      out << "  Eta range             : [" << eta_range[0] << ", "
          << eta_range[1] << "]\n";
    } else {
      const auto& theta_range = cfg.theta_range();
      out << "  Theta range           : ["
          << theta_range[0] / detray::unit<scalar_t>::rad << ", "
          << theta_range[1] / detray::unit<scalar_t>::rad << ") [rad]\n";
    }

    // Origin
    out << "  Origin                : [" << ori[0] / detray::unit<scalar_t>::mm
        << ", " << ori[1] / detray::unit<scalar_t>::mm << ", "
        << ori[2] / detray::unit<scalar_t>::mm << "] [mm]\n";

    return out;
  }
};

}  // namespace detray
