// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// Detray test include(s)
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"
#include "detray/test/framework/test_configuration.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <limits>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Configuration for the material ray scan.
///
/// @tparam algebra_t the algebra type of the detector that is scanned
template <concepts::algebra algebra_t>
struct material_scan_config
    : public detray::test::configuration<dscalar<algebra_t>> {
  /// Scalar type used by the scan.
  using scalar_type = dscalar<algebra_t>;
  /// Configuration type of the ray generator.
  using trk_gen_config_t = uniform_track_generator_config<scalar_type>;

  std::string m_name{"material_scan"};
  trk_gen_config_t m_trk_gen_cfg{};
  /// Name of the output file, containing the complete ray material traces
  std::string m_material_file{"material_scan"};
  /// Perform overlaps removal (needed for detector converted from ACTS)
  bool m_overlaps_removal{true};
  /// Tolerance for overlaps
  float m_overlaps_tol{1e-4f * unit<float>::mm};
  /// Read the truth material from the material file instead of running the
  /// ray scan
  bool m_truth_from_file{false};

  /// Getters
  /// @{
  const std::string &name() const { return m_name; }
  trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
  const trk_gen_config_t &track_generator() const { return m_trk_gen_cfg; }
  const std::string &material_file() const { return m_material_file; }
  bool overlaps_removal() const { return m_overlaps_removal; }
  float overlaps_tol() const { return m_overlaps_tol; }
  /// @returns whether the truth material is read from file
  bool truth_from_file() const { return m_truth_from_file; }
  /// @}

  /// Setters
  /// @{
  material_scan_config &name(const std::string &n) {
    m_name = n;
    return *this;
  }
  material_scan_config &material_file(const std::string &f) {
    m_material_file = f;
    return *this;
  }
  material_scan_config &overlaps_removal(const bool o) {
    m_overlaps_removal = o;
    return *this;
  }
  material_scan_config &overlaps_tol(const scalar_type tol) {
    m_overlaps_tol = static_cast<float>(tol);
    return *this;
  }
  /// Read the truth material from file instead of running the ray scan
  /// @param b read from file if true
  /// @returns this configuration
  material_scan_config &truth_from_file(const bool b) {
    m_truth_from_file = b;
    return *this;
  }
  /// @}
};

/// @brief Configuration for a detector scan test.
template <concepts::algebra algebra_t>
struct material_validation_config
    : public detray::test::configuration<dscalar<algebra_t>> {
  /// Scalar type used by the validation.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional vector type used by the validation.
  using vector3_type = dvector3D<algebra_t>;
  /// Base test configuration type.
  using base_type = detray::test::configuration<scalar_type>;

  /// Name of the test
  std::string m_name{"material_validation"};
  /// Vecmem memory resource for the device allocations
  vecmem::memory_resource *m_dev_mr{nullptr};
  /// Name of the output file, containing the complete ray material traces
  std::string m_material_file{"navigation_material_trace"};
  /// The maximal number of test tracks to run
  std::size_t m_n_tracks{detray::detail::invalid_value<std::size_t>()};
  /// Allowed relative discrepancy between truth and navigation material
  scalar_type m_rel_error{0.001f};
  /// Additional material steps per track on top of the steps found by the
  /// material scan (capacity of the device material step buffers)
  std::size_t m_step_contingency{10u};

  /// Getters
  /// @{
  /// Get the test name.
  /// @return The configured test name.
  const std::string &name() const { return m_name; }
  /// Get the device memory resource.
  /// @return The device memory resource.
  vecmem::memory_resource *device_mr() const { return m_dev_mr; }
  /// Get the material trace file name.
  /// @return The material trace file name.
  const std::string &material_file() const { return m_material_file; }
  /// Get the maximum number of test tracks.
  /// @return The maximum number of test tracks.
  std::size_t n_tracks() const { return m_n_tracks; }
  /// Get the allowed relative discrepancy.
  /// @return The allowed relative discrepancy.
  scalar_type relative_error() const { return m_rel_error; }
  /// Get the additional material steps per track.
  /// @return The number of additional material steps.
  std::size_t step_contingency() const { return m_step_contingency; }
  /// @}

  /// Setters
  /// @{
  /// Set the test name.
  /// @param n New test name.
  /// @return This configuration.
  material_validation_config &name(const std::string &n) {
    m_name = n;
    return *this;
  }
  /// Set the device memory resource.
  /// @param mr Memory resource for device allocations.
  /// @return This configuration.
  material_validation_config &device_mr(vecmem::memory_resource *mr) {
    m_dev_mr = mr;
    return *this;
  }
  /// Set the material trace file name.
  /// @param f New file name.
  /// @return This configuration.
  material_validation_config &material_file(const std::string &f) {
    m_material_file = f;
    return *this;
  }
  /// Set the maximum number of test tracks.
  /// @param n Number of tracks.
  /// @return This configuration.
  material_validation_config &n_tracks(std::size_t n) {
    m_n_tracks = n;
    return *this;
  }
  /// Set the allowed relative discrepancy.
  /// @param re Relative discrepancy.
  /// @return This configuration.
  material_validation_config &relative_error(scalar_type re) {
    assert(re > 0.f);
    m_rel_error = re;
    return *this;
  }
  /// Set the additional material steps per track.
  /// @param n Number of additional material steps.
  /// @return This configuration.
  material_validation_config &step_contingency(std::size_t n) {
    m_step_contingency = n;
    return *this;
  }
  /// @}
};

}  // namespace detray::test
