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
#include "detray/geometry/surface.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/tools/event_generator/uniform_track_generator.hpp"
#include "detray/test/utils/data_record.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/material_recording.hpp"

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

/// @brief Configuration for the material ray scan.
///
/// @tparam algebra_t the algebra type of the detector that is scanned
template <concepts::algebra algebra_t>
struct material_scan_config
    : public detray::test::configuration<dscalar<algebra_t>,
                                         material_scan_config<algebra_t>> {
  /// Scalar type used by the scan.
  using scalar_type = dscalar<algebra_t>;
  /// Configuration type of the ray generator.
  using trk_gen_config_t = uniform_track_generator_config<scalar_type>;
  /// Base test configuration type.
  using base_type =
      detray::test::configuration<scalar_type, material_scan_config>;

  /// Default constructor calling base constructor with specific name
  material_scan_config() : base_type("material_scan") {};

  /// The configuration of the uniform track generator
  trk_gen_config_t m_trk_gen_cfg{};
  /// Name of the output file, containing the complete ray material traces
  std::string m_material_file{"material_scan"};
  /// Tolerance for overlaps
  float m_overlaps_tol{1e-4f * unit<float>::mm};
  /// Perform overlaps removal (needed for detector converted from ACTS)
  bool m_overlaps_removal{true};
  /// Read the scan data from the material file instead of running the ray scan
  bool m_read_from_file{false};

  /// Getters
  /// @{
  /// Get the track generator configuration.
  /// @return the track generator configuration.
  trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
  /// Get the track generator configuration.
  /// @return the const track generator configuration.
  const trk_gen_config_t &track_generator() const { return m_trk_gen_cfg; }
  /// Get the material trace file name.
  /// @return The material trace file name.
  const std::string &material_file() const { return m_material_file; }
  /// Get the allowed overlaps tolerance.
  /// @return The allowed relative discrepancy.
  float overlaps_tol() const { return m_overlaps_tol; }
  /// Whether to perform overlaps removal during the scan.
  /// @return The overlaps removal toggle.
  bool overlaps_removal() const { return m_overlaps_removal; }
  /// @returns whether the truth material is read from file
  bool read_from_file() const { return m_read_from_file; }
  /// @}

  /// Setters
  /// @{
  /// Set the material trace file name.
  /// @param f New file name.
  /// @return This configuration.
  material_scan_config &material_file(const std::string &f) {
    m_material_file = f;
    return *this;
  }
  /// Set the allowed overlaps tolerance.
  /// @param tol The allowed relative discrepancy.
  /// @return This configuration.
  material_scan_config &overlaps_tol(const scalar_type tol) {
    m_overlaps_tol = static_cast<float>(tol);
    return *this;
  }
  /// Toggle whether to perform overlaps removal during the scan.
  /// @param o The overlaps removal toggle.
  /// @return This configuration.
  material_scan_config &overlaps_removal(const bool o) {
    m_overlaps_removal = o;
    return *this;
  }
  /// Read the truth material from file instead of running the ray scan
  /// @param b read from file if true
  /// @returns this configuration
  material_scan_config &read_from_file(const bool b) {
    m_read_from_file = b;
    return *this;
  }
  /// @}
};

/// Read the material scan data from file and generate the corresponding tracks
///
/// @param file_name the name of the file that contains the scan data
/// @param cfg the material scan configuration that was used to write the file
///
/// @returns the truth tracks and material
template <concepts::algebra algebra_t>
inline auto read_material_scan(const std::string &file_name,
                               const material_scan_config<algebra_t> &cfg)
    -> material_trace_data<algebra_t> {
  using scalar_t = dscalar<algebra_t>;

  // Fill trace data from file
  material_trace_data<algebra_t> scan_data{};
  scan_data.track_materials = read_track_material(file_name);
  const dvector<track_material> &track_materials = scan_data.track_materials;

  // Reconstruct the initial free track parameters for the scan ray
  scan_data.tracks.reserve(track_materials.size());
  for (const track_material &trk_mat : track_materials) {
    const double phi{trk_mat.phi};
    const double theta{2. * math::atan(math::exp(-trk_mat.eta))};
    const double sin_theta{math::sin(theta)};

    // No momentum needed for rays
    const dvector3D<algebra_t> dir{
        static_cast<scalar_t>(math::cos(phi) * sin_theta),
        static_cast<scalar_t>(math::sin(phi) * sin_theta),
        static_cast<scalar_t>(math::cos(theta))};

    scan_data.tracks.push_back(
        {cfg.track_generator().origin(), 0.f, vector::normalize(dir), 0.f});
  }

  std::clog << "-----------------------------------\n"
            << "Read " << scan_data.tracks.size() << " tracks from "
            << file_name << "\n"
            << "-----------------------------------\n"
            << std::endl;

  return scan_data;
}

namespace material_scanner {

/// Run the material ray scan on a given detector and write the truth material
/// to file (or read it from file, if configured)
///
/// @param gctx the geometry context to scan
/// @param det the detector to scan
/// @param names the volume names of the detector
/// @param cfg the material scan configuration
///
/// @returns the truth tracks and material
template <concepts::detector detector_t>
inline auto run(
    const typename detector_t::geometry_context &gctx, const detector_t &det,
    const typename detector_t::name_map &names,
    const material_scan_config<typename detector_t::algebra_type> &cfg)
    -> material_trace_data<typename detector_t::algebra_type> {
  using algebra_t = typename detector_t::algebra_type;
  using point2_t = dpoint2D<algebra_t>;
  using ray_t = detail::ray<algebra_t>;
  using track_generator_t = uniform_track_generator<ray_t>;

  const std::string det_name{det.name(names)};
  const std::string file_name =
      material_file_path(det_name, cfg.material_file()).string();

  if (cfg.read_from_file()) {
    return read_material_scan(file_name, cfg);
  }

  std::size_t n_tracks{0u};
  auto ray_generator = track_generator_t(cfg.track_generator());

  DETRAY_INFO_HOST("Running material scan on: " << det_name << "\n("
                                                << ray_generator.size()
                                                << " rays) ...\n");

  // Trace material per ray
  material_trace_data<algebra_t> scan_data{};
  scan_data.tracks.reserve(ray_generator.size());
  scan_data.track_materials.reserve(ray_generator.size());
  scan_data.material_steps.reserve(ray_generator.size());

  for (const auto ray : ray_generator) {
    // Record all intersections and surfaces along the ray
    auto intersection_record =
        detector_scanner::run<detray::ray_scan>(gctx, det, ray);

    // Remove certain allowed duplications
    if (cfg.overlaps_removal()) {
      detector_scanner::overlaps_removal(intersection_record,
                                         cfg.overlaps_tol());
    }

    if (intersection_record.empty()) {
      DETRAY_FATAL_HOST("Intersection trace empty for ray "
                        << n_tracks << "/" << ray_generator.size() << ": "
                        << ray);
      break;
    }

    // Record track parameters
    scan_data.tracks.push_back({ray.pos(), 0.f, ray.dir(), 0.f});

    // New track material entry
    track_material track_mat{};
    track_mat.eta = vector::eta(ray.dir());
    track_mat.phi = vector::phi(ray.dir());

    // New material records collection (matertial steps)
    dvector<material_record> mat_steps{};
    if (!scan_data.material_steps.empty()) {
      mat_steps.reserve(scan_data.material_steps.back().size());
    }

    // Record material for this ray
    for (const auto &[i, record] :
         detray::views::enumerate(intersection_record)) {
      // Check whether this record has a successor
      if (i < intersection_record.size() - 1) {
        const auto &current_intr = record.intersection();
        const auto &next_intr = intersection_record[i + 1].intersection();

        const bool is_same_intrs{next_intr == current_intr};
        const bool current_is_pt{current_intr.surface().is_portal()};
        const bool next_is_pt{next_intr.surface().is_portal()};

        const bool is_dummy_record{i == 0};

        // Prevent double counting of material on adjacent portals
        // (the navigator automatically skips the exit portal)
        if ((is_same_intrs && next_is_pt && current_is_pt) || is_dummy_record) {
          continue;
        }
      }

      // Don't count the last portal, because navigation terminates
      // before the material is counted
      if (detail::is_invalid_value(record.intersection().volume_link())) {
        continue;
      }

      const auto sf = geometry::surface{det, record.intersection().surface()};

      if (!sf.has_material()) {
        continue;
      }

      const auto &p = record.intersection().local();
      const auto mat_record = sf.template visit_material<get_material_record>(
          point2_t{p[0], p[1]}, cos_angle(gctx, sf, ray.dir(), p));

      const double seg{mat_record.path};
      const double t{mat_record.thickness};
      const double mx0{mat_record.mat_X0};
      const double ml0{mat_record.mat_L0};

      if (mx0 > 0.f) {
        track_mat.sX0 += seg / mx0;
        track_mat.tX0 += t / mx0;
      } else {
        DETRAY_ERROR_HOST(
            "Encountered invalid X_0: " << mx0 << "\nOn surface: " << sf);
        continue;
      }
      if (ml0 > 0.f) {
        track_mat.sL0 += seg / ml0;
        track_mat.tL0 += t / ml0;
      } else {
        DETRAY_ERROR_HOST(
            "Encountered invalid L_0: " << ml0 << "\nOn surface: " << sf);
        continue;
      }

      // Record the material step for this track
      mat_steps.push_back(mat_record);
    }

    if (track_mat.sX0 == 0.f || track_mat.sL0 == 0.f || track_mat.tX0 == 0.f ||
        track_mat.tL0 == 0.f) {
      DETRAY_VERBOSE_HOST("No material recorded for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray);
    }

    scan_data.track_materials.push_back(track_mat);
    scan_data.material_steps.push_back(std::move(mat_steps));

    ++n_tracks;
  }

  std::clog << "-----------------------------------\n"
            << "Tested " << n_tracks << " tracks\n"
            << "-----------------------------------\n"
            << std::endl;

  // Write recorded material to csv file
  write_track_material(file_name, scan_data.track_materials);

  return scan_data;
}

}  // namespace material_scanner

}  // namespace detray::validation
