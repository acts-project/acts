// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/utils/data_record.hpp"
#include "detray/test/validation/material_comparison.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace detray::validation {

/// @brief Configuration for a detector scan test.
template <concepts::algebra algebra_t>
struct material_validation_config
    : public detray::test::configuration<
          dscalar<algebra_t>, material_validation_config<algebra_t>> {
  /// Scalar type used by the validation.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional vector type used by the validation.
  using vector3_type = dvector3D<algebra_t>;
  /// Base test configuration type.
  using base_type =
      detray::test::configuration<scalar_type, material_validation_config>;

  /// Default constructor calling base constructor with specific name
  material_validation_config() : base_type("material_validation") {};

  /// Vecmem memory resource for the device allocations
  vecmem::memory_resource *m_dev_mr{nullptr};
  /// Name of the output file, containing the complete ray material traces
  std::string m_material_file{"navigation_material_trace"};
  /// Allowed relative discrepancy between truth and navigation material
  scalar_type m_rel_error{0.001f};
  /// The maximal number of test tracks to run
  std::size_t m_n_tracks{detray::detail::invalid_value<std::size_t>()};
  /// Additional material steps per track on top of the steps found by the
  /// material scan (capacity of the device material step buffers)
  std::size_t m_step_contingency{10u};

  /// Getters
  /// @{
  /// Get the device memory resource.
  /// @return The device memory resource.
  vecmem::memory_resource *device_mr() const { return m_dev_mr; }
  /// Get the material trace file name.
  /// @return The material trace file name.
  const std::string &material_file() const { return m_material_file; }
  /// Get the allowed relative discrepancy.
  /// @return The allowed relative discrepancy.
  scalar_type relative_error() const { return m_rel_error; }
  /// Get the maximum number of test tracks.
  /// @return The maximum number of test tracks.
  std::size_t n_tracks() const { return m_n_tracks; }
  /// Get the additional material steps per track.
  /// @return The number of additional material steps.
  std::size_t step_contingency() const { return m_step_contingency; }
  /// @}

  /// Setters
  /// @{
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
  /// Set the allowed relative discrepancy.
  /// @param re Relative discrepancy.
  /// @return This configuration.
  material_validation_config &relative_error(scalar_type re) {
    assert(re > 0.f);
    m_rel_error = re;
    return *this;
  }
  /// Set the maximum number of test tracks.
  /// @param n Number of tracks.
  /// @return This configuration.
  material_validation_config &n_tracks(std::size_t n) {
    m_n_tracks = n;
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

namespace material_validator {

/// @brief Result of the material validation
struct result {
  /// Number of tracks that were compared
  std::size_t n_tracks{0u};
  /// Number of tracks with a difference in integrated material
  std::size_t n_diff_total_mat{0u};
  /// Number of tracks with a difference in material steps
  std::size_t n_diff_mat_steps{0u};

  /// @returns true if no material discrepancy was found
  bool passed() const {
    return (n_diff_total_mat == 0u && n_diff_mat_steps == 0u);
  }
};

/// Trace the material during propagation and compare it to the truth material
///
/// Prints an error message for every material discrepancy and writes the
/// recorded material to file.
///
/// @tparam backend_t functor that runs the propagation and records the
///                   material (e.g. on host or device)
///
/// @param host_mr the host memory resource
/// @param det the detector to be validated
/// @param names the volume names of the detector
/// @param truth_traces the truth tracks and their collected material
/// @param cfg the material comparison configuration
///
/// @returns the number of compared and failed tracks
template <typename material_recorder_t, concepts::detector detector_t>
inline result run(
    vecmem::memory_resource *host_mr, const detector_t &det,
    const typename detector_t::name_map &names,
    const material_trace_data<typename detector_t::algebra_type> &truth_traces,
    const material_validation_config<typename detector_t::algebra_type> &cfg) {
  using algebra_t = typename detector_t::algebra_type;

  const auto &tracks = truth_traces.tracks;
  const auto &truth_track_mat = truth_traces.track_materials;
  const auto &truth_mat_steps = truth_traces.material_steps;

  assert(tracks.size() == truth_track_mat.size());
  assert(tracks.size() == truth_mat_steps.size());

  const std::string det_name{det.name(names)};

  DETRAY_INFO_HOST("Running material validation on: " << det_name << "...\n");

  // Capacities of the device material step allocations: Number of material
  // steps found by the material scan plus contingency
  constexpr std::size_t default_capacity{80u};
  std::vector<std::size_t> capacities(tracks.size(), default_capacity);
  if (truth_traces.material_steps.size() == tracks.size()) {
    for (std::size_t trk_idx = 0u; trk_idx < tracks.size(); ++trk_idx) {
      capacities.at(trk_idx) = truth_traces.material_steps.at(trk_idx).size() +
                               cfg.step_contingency();
    }
  }

  // Run the propagation and record the material (accumulated and steps)
  constexpr material_recorder_t material_recorder{};
  material_trace_data<algebra_t> recorded_traces = material_recorder(
      host_mr, cfg.device_mr(), det, cfg.propagation(), tracks, capacities);

  const auto &rec_tracks = recorded_traces.tracks;
  const auto &rec_track_mat = recorded_traces.track_materials;
  const auto &rec_mat_steps = recorded_traces.material_steps;

  // Check that the material recording did not run into any errors
  if (tracks.size() != rec_tracks.size()) {
    throw std::runtime_error(
        "Material comparison: Number of material records (" +
        std::to_string(rec_track_mat.size()) +
        ") does not match the number of tracks (" +
        std::to_string(tracks.size()) + ")");
  }
  assert(rec_tracks.size() == rec_track_mat.size());
  assert(rec_tracks.size() == rec_mat_steps.size());

  // Compare the traces
  result comp_result{};
  for (std::size_t i = 0u; i < tracks.size(); ++i) {
    // Only run comparison for the configured number of tracks
    if (comp_result.n_tracks >= cfg.n_tracks()) {
      break;
    }

    // Compare the single material steps and total integrated material per track
    constexpr bool verbose{true};
    const auto [bad_steps, bad_total] = validation::compare_material_traces(
        truth_mat_steps.at(i), truth_track_mat.at(i), rec_mat_steps.at(i),
        rec_track_mat.at(i), i, cfg.relative_error(), verbose);

    // Total accumulated material differs (e.g. different number of mat. steps)
    if (bad_total) {
      ++comp_result.n_diff_total_mat;
    }
    // Different material seen at at least one material step
    if (bad_steps) {
      ++comp_result.n_diff_mat_steps;
    }

    ++comp_result.n_tracks;
  }

  std::clog << "-----------------------------------\n"
            << "Tested " << comp_result.n_tracks << " tracks\n"
            << "-----------------------------------\n"
            << std::endl;

  // Write accumulated material per track to file
  write_track_material(material_file_path(det_name, cfg.material_file(),
                                          material_recorder_t::name)
                           .string(),
                       rec_track_mat);

  return comp_result;
}

}  // namespace material_validator

}  // namespace detray::validation
