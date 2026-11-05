// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/pdg_particle.hpp"

// detray plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"

// Detray test include(s)
#include "detray/test/framework/test_configuration.hpp"

// System include(s)
#include <limits>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Configuration for a detector scan test.
template <concepts::algebra algebra_t>
struct navigation_validation_config
    : public detray::test::configuration<dscalar<algebra_t>> {
  /// Scalar type used by the validation.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional vector type used by the validation.
  using vector3_type = dvector3D<algebra_t>;
  /// Base test configuration type.
  using base_type = detray::test::configuration<scalar_type>;

  /// Name of the test
  std::string m_name{"navigation_validation"};
  /// Name of the input file, containing the complete ray scan traces
  std::string m_intersection_file{"truth_intersections.csv"};
  /// Name of the input file containing track parameters.
  std::string m_track_param_file{"truth_trk_parameters.csv"};
  /// The maximum number of test tracks to run
  std::size_t m_n_tracks{detray::detail::invalid_value<std::size_t>()};
  /// Particle hypothesis (truth particle from simulation)
  pdg_particle<scalar_type> m_ptc_hypo{muon<scalar_type>()};
  /// Navigation direction
  navigation::direction m_nav_dir{navigation::direction::e_forward};
  /// Collect only the sensitive intersections for comparison
  bool m_collect_sensitives_only{false};
  /// Toggle SVG display
  bool m_display_svg{false};
  /// Drop an SVG only if navigation missed a surface
  bool m_display_only_missed{false};
  /// Whether to stop execution at the first error
  bool m_fail_on_diff{true};
  /// Verbosity of the console output
  bool m_verbose{true};
  /// Configured momentum range of the test sample (only needed to generate
  /// correct file names). If none was passed, it will be determined from
  /// the track data (imprecise!)
  darray<scalar_type, 2> m_p_range{
      detray::detail::invalid_value<scalar_type>(),
      detray::detail::invalid_value<scalar_type>()};
  /// B-field vector for helix
  vector3_type m_B{0.f * unit<scalar_type>::T, 0.f * unit<scalar_type>::T,
                   2.f * unit<scalar_type>::T};
  /// Visualization style to be applied to the SVGs
  detray::svgtools::styling::style m_style =
      detray::svgtools::styling::tableau_colorblind::style;

  /// Getters
  /// @{
  /// Get the test name.
  /// @return The configured test name.
  const std::string &name() const { return m_name; }
  /// Get the intersection input file name.
  /// @return The intersection input file name.
  const std::string &intersection_file() const { return m_intersection_file; }
  /// Get the track parameter input file name.
  /// @return The track parameter input file name.
  const std::string &track_param_file() const { return m_track_param_file; }
  /// Get the maximum number of tracks.
  /// @return The maximum number of tracks.
  std::size_t n_tracks() const { return m_n_tracks; }
  /// Get the particle hypothesis.
  /// @return The particle hypothesis.
  pdg_particle<scalar_type> ptc_hypothesis() const { return m_ptc_hypo; }
  /// Get the navigation direction.
  /// @return The navigation direction.
  navigation::direction navigation_direction() const { return m_nav_dir; }
  /// Check whether only sensitive surfaces are collected.
  /// @return Whether only sensitive surfaces are collected.
  bool collect_sensitives_only() const { return m_collect_sensitives_only; }
  /// Check whether SVG output is enabled.
  /// @return Whether SVG output is enabled.
  bool display_svg() const { return m_display_svg; }
  /// Check whether SVG output is limited to missed surfaces.
  /// @return Whether SVG output is limited to missed surfaces.
  bool display_only_missed() const { return m_display_only_missed; }
  /// Check whether validation stops at the first difference.
  /// @return Whether validation stops at the first difference.
  bool fail_on_diff() const { return m_fail_on_diff; }
  /// Check whether verbose output is enabled.
  /// @return Whether verbose output is enabled.
  bool verbose() const { return m_verbose; }
  /// Get the configured momentum range.
  /// @return The configured momentum range.
  darray<scalar_type, 2> p_range() const { return m_p_range; }
  /// Get the magnetic field vector.
  /// @return The magnetic field vector.
  const vector3_type &B_vector() const { return m_B; }
  /// Get the SVG style.
  /// @return The SVG style.
  const auto &svg_style() const { return m_style; }
  /// @}

  /// Setters
  /// @{
  /// Set the test name.
  /// @param n New test name.
  /// @return This configuration.
  navigation_validation_config &name(const std::string &n) {
    m_name = n;
    return *this;
  }
  /// Set the intersection input file name.
  /// @param f New file name.
  /// @return This configuration.
  navigation_validation_config &intersection_file(const std::string &f) {
    m_intersection_file = f;
    return *this;
  }
  /// Set the track parameter input file name.
  /// @param f New file name.
  /// @return This configuration.
  navigation_validation_config &track_param_file(const std::string &f) {
    m_track_param_file = f;
    return *this;
  }
  /// Set the particle hypothesis.
  /// @param pdg_ptc New particle hypothesis.
  /// @return This configuration.
  navigation_validation_config &ptc_hypothesis(
      pdg_particle<scalar_type> pdg_ptc) {
    m_ptc_hypo = pdg_ptc;
    return *this;
  }
  /// Set the navigation direction.
  /// @param dir New direction.
  /// @return This configuration.
  navigation_validation_config &navigation_direction(
      const navigation::direction dir) {
    m_nav_dir = dir;
    return *this;
  }
  /// Set whether to collect only sensitive surfaces.
  /// @param only_sensitives Whether to restrict collection.
  /// @return This configuration.
  navigation_validation_config &collect_sensitives_only(
      const bool only_sensitives) {
    m_collect_sensitives_only = only_sensitives;
    return *this;
  }
  /// Set whether SVG output is enabled.
  /// @param do_display Whether to display SVGs.
  /// @return This configuration.
  navigation_validation_config &display_svg(const bool do_display) {
    m_display_svg = do_display;
    return *this;
  }
  /// Set whether SVG output is limited to missed surfaces.
  /// @param only_missed Whether to display only misses.
  /// @return This configuration.
  navigation_validation_config &display_only_missed(const bool only_missed) {
    m_display_only_missed = only_missed;
    return *this;
  }
  /// Set whether validation stops at the first difference.
  /// @param fail_on_diff Whether to stop on a difference.
  /// @return This configuration.
  navigation_validation_config &fail_on_diff(const bool fail_on_diff) {
    m_fail_on_diff = fail_on_diff;
    return *this;
  }
  /// Set whether verbose output is enabled.
  /// @param v Whether to enable verbose output.
  /// @return This configuration.
  navigation_validation_config &verbose(const bool v) {
    m_verbose = v;
    return *this;
  }
  /// Set the maximum number of tracks.
  /// @param n Number of tracks.
  /// @return This configuration.
  navigation_validation_config &n_tracks(std::size_t n) {
    m_n_tracks = n;
    return *this;
  }
  /// Set the momentum range.
  /// @param pr New momentum range.
  /// @return This configuration.
  navigation_validation_config &p_range(const darray<scalar_type, 2> pr) {
    m_p_range = pr;
    return *this;
  }
  /// Set the magnetic field vector.
  /// @param B New field vector.
  /// @return This configuration.
  navigation_validation_config &B_vector(const vector3_type &B) {
    m_B = B;
    return *this;
  }
  /// Set the magnetic field vector components.
  /// @param B0 First field component.
  /// @param B1 Second field component.
  /// @param B2 Third field component.
  /// @return This configuration.
  navigation_validation_config &B_vector(const scalar_type B0,
                                         const scalar_type B1,
                                         const scalar_type B2) {
    m_B = vector3_type{B0, B1, B2};
    return *this;
  }
  /// @}
};

}  // namespace detray::test
