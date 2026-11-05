// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"

// Detray test include(s)
#include "detray/test/framework/test_configuration.hpp"

// System include(s)
#include <limits>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Configuration for a detector scan test.
template <typename track_generator_t, concepts::algebra algebra_t>
struct detector_scan_config
    : public detray::test::configuration<dscalar<algebra_t>> {
  /// Scalar type used by the scan.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional vector type used by the scan.
  using vector3_type = dvector3D<algebra_t>;
  /// Base test configuration type.
  using base_type = detray::test::configuration<scalar_type>;
  /// Track generator configuration type.
  using trk_gen_config_t = typename track_generator_t::configuration;

  /// Name of the test
  std::string m_name{"detector_scan"};
  /// Name of the input file, containing the complete ray scan traces
  std::string m_intersection_file{"truth_intersections"};
  /// Name of the input file containing track parameters.
  std::string m_track_param_file{"truth_trk_parameters"};
  /// Mask tolerance for the intersectors
  scalar_type m_mask_tol{std::numeric_limits<scalar_type>::epsilon()};
  /// B-field vector for helix
  vector3_type m_B{0.f * unit<scalar_type>::T, 0.f * unit<scalar_type>::T,
                   2.f * unit<scalar_type>::T};
  /// Configuration of the ray generator
  trk_gen_config_t m_trk_gen_cfg{};
  /// Perform overlaps removal (needed for detector converted from ACTS)
  bool m_overlaps_removal{true};
  /// Tolerance for overlaps
  float m_overlaps_tol{1e-4f * unit<float>::mm};
  /// Write intersection points for plotting
  bool m_write_inters{false};
  /// Visualization style to be applied to the svgs
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
  /// Get the intersector mask tolerance.
  /// @return The intersector mask tolerance.
  scalar_type mask_tolerance() const { return m_mask_tol; }
  /// Get the magnetic field vector.
  /// @return The magnetic field vector.
  const vector3_type &B_vector() { return m_B; }
  /// Check whether intersections are written.
  /// @return Whether intersections are written.
  bool write_intersections() const { return m_write_inters; }
  /// Get the mutable track generator configuration.
  /// @return The mutable track generator configuration.
  trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
  /// Get the track generator configuration.
  /// @return The track generator configuration.
  const trk_gen_config_t &track_generator() const { return m_trk_gen_cfg; }
  /// Check whether overlap removal is enabled.
  /// @return Whether overlap removal is enabled.
  bool overlaps_removal() const { return m_overlaps_removal; }
  /// Get the overlap tolerance.
  /// @return The overlap tolerance.
  float overlaps_tol() const { return m_overlaps_tol; }
  /// Get the SVG style.
  /// @return The SVG style.
  const auto &svg_style() const { return m_style; }
  /// @}

  /// Setters
  /// @{
  /// Set the test name.
  /// @param n New test name.
  /// @return This configuration.
  detector_scan_config &name(const std::string_view n) {
    m_name = n;
    return *this;
  }
  /// Set the intersection input file name.
  /// @param f New file name.
  /// @return This configuration.
  detector_scan_config &intersection_file(const std::string_view f) {
    m_intersection_file = f;
    return *this;
  }
  /// Set the track parameter input file name.
  /// @param f New file name.
  /// @return This configuration.
  detector_scan_config &track_param_file(const std::string_view f) {
    m_track_param_file = f;
    return *this;
  }
  /// Set the intersector mask tolerance.
  /// @param tol New tolerance.
  /// @return This configuration.
  detector_scan_config &mask_tolerance(const scalar_type tol) {
    m_mask_tol = tol;
    return *this;
  }
  /// Set the magnetic field vector.
  /// @param B New field vector.
  /// @return This configuration.
  detector_scan_config &B_vector(const vector3_type &B) {
    m_B = B;
    return *this;
  }
  /// Set the magnetic field vector components.
  /// @param B0 First field component.
  /// @param B1 Second field component.
  /// @param B2 Third field component.
  /// @return This configuration.
  detector_scan_config &B_vector(const scalar_type B0, const scalar_type B1,
                                 const scalar_type B2) {
    m_B = vector3_type{B0, B1, B2};
    return *this;
  }
  /// Set whether overlap removal is enabled.
  /// @param o Whether to remove overlaps.
  /// @return This configuration.
  detector_scan_config &overlaps_removal(const bool o) {
    m_overlaps_removal = o;
    return *this;
  }
  /// Set the overlap tolerance.
  /// @param tol New tolerance.
  /// @return This configuration.
  detector_scan_config &overlaps_tol(const scalar_type tol) {
    m_overlaps_tol = static_cast<float>(tol);
    return *this;
  }
  /// Set whether intersections are written.
  /// @param do_write Whether to write intersections.
  /// @return This configuration.
  detector_scan_config &write_intersections(const bool do_write) {
    m_write_inters = do_write;
    return *this;
  }
  /// @}
};

}  // namespace detray::test
