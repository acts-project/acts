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
#include "detray/test/common/track_generators.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/material_validation_config.hpp"
#include "detray/test/validation/material_validation_utils.hpp"

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

namespace detray::test {

/// Run the material validation on the host
struct run_material_validation {
  static constexpr std::string_view name{"cpu"};

  template <concepts::detector detector_t>
  auto operator()(
      vecmem::memory_resource *host_mr, vecmem::memory_resource * /*mr*/,
      const detector_t &det, const propagation::config &cfg,
      const dvector<free_track_parameters<typename detector_t::algebra_type>>
          &tracks,
      const std::vector<std::size_t> & /*mask*/ = {}) {
    typename detector_t::geometry_context gctx{};

    dvector<material_validator::track_material> track_mat_vec{host_mr};
    track_mat_vec.reserve(tracks.size());

    dvector<dvector<material_record>> mat_steps_vec{host_mr};
    mat_steps_vec.reserve(tracks.size());

    for (const auto &[i, track] : detray::views::enumerate(tracks)) {
      auto [success, track_mat, mat_steps] =
          detray::material_validator::record_material(gctx, host_mr, det, cfg,
                                                      track);
      track_mat_vec.push_back(track_mat);
      mat_steps_vec.push_back(std::move(mat_steps));

      if (!success) {
        DETRAY_ERROR_HOST("Propagation failed for track "
                          << i << ": "
                          << "Material record may be incomplete!");
      }
    }

    return std::make_tuple(std::move(track_mat_vec), std::move(mat_steps_vec));
  }
};

}  // namespace detray::test

namespace detray::material_validator {

/// @brief Truth data of the material scan
///
/// @tparam algebra_t the algebra type of the scanned detector
template <concepts::algebra algebra_t>
struct material_scan_data {
  /// Test tracks (one per scanned ray)
  dvector<free_track_parameters<algebra_t>> tracks{};
  /// Accumulated truth material per track
  dvector<track_material> material{};
  /// Number of surfaces with material per track (empty if the truth data was
  /// read from file)
  std::vector<std::size_t> n_steps{};
};

/// @brief Result of the material validation
struct material_validation_result {
  /// Number of tracks that were compared
  std::size_t n_tracks{0u};
  /// Number of tracks with a material discrepancy
  std::size_t n_failed{0u};

  /// @returns true if no material discrepancy was found
  bool passed() const { return n_failed == 0u; }
};

/// Name of a material data file: "<detector>_<file stem>[_<suffix>].csv"
///
/// @param det_name the name of the detector
/// @param material_file the configured output file (only the parent path and
///                      the stem are used)
/// @param suffix optional suffix, e.g. the name of the compute backend
///
/// @returns the full path to the material data file
inline std::filesystem::path material_file_path(const std::string &det_name,
                                                const std::string &material_file,
                                                std::string_view suffix = {}) {
  const std::filesystem::path mat_path{material_file};

  std::string file_name{det_name + "_" + mat_path.stem().string()};
  if (!suffix.empty()) {
    file_name += "_" + std::string(suffix);
  }

  return mat_path.parent_path() / (file_name + ".csv");
}

/// Read the truth material from file and generate the corresponding tracks
///
/// The ray generator is not costly to run, so the tracks are generated again
/// and checked against the eta and phi values in the material file.
///
/// @param det_name the name of the detector
/// @param cfg the material scan configuration that was used to write the file
///
/// @returns the truth tracks and material
template <concepts::algebra algebra_t>
inline auto read_material_scan(const std::string &det_name,
                               const test::material_scan_config<algebra_t> &cfg)
    -> material_scan_data<algebra_t> {
  using ray_t = detail::ray<algebra_t>;
  using track_generator_t = uniform_track_generator<ray_t>;

  // Tolerance for the eta and phi values in the csv file
  constexpr double tol{1e-4};

  const auto file_name{material_file_path(det_name, cfg.material_file())};

  material_scan_data<algebra_t> scan_data{};
  scan_data.material = read_material(file_name.string());

  auto ray_generator = track_generator_t(cfg.track_generator());
  if (scan_data.material.size() > ray_generator.size()) {
    throw std::invalid_argument(
        "Material file " + file_name.string() + " contains more tracks (" +
        std::to_string(scan_data.material.size()) +
        ") than the material scan configuration (" +
        std::to_string(ray_generator.size()) + ")");
  }

  scan_data.tracks.reserve(scan_data.material.size());
  for (const auto ray : ray_generator) {
    const std::size_t i{scan_data.tracks.size()};
    if (i == scan_data.material.size()) {
      break;
    }

    const auto &rec = scan_data.material[i];
    const double eta{vector::eta(ray.dir())};
    const double phi{vector::phi(ray.dir())};
    if (math::fabs(rec.eta - eta) > tol || math::fabs(rec.phi - phi) > tol) {
      throw std::invalid_argument(
          "Material file " + file_name.string() +
          " does not match the material scan configuration: Track " +
          std::to_string(i) + " has (eta, phi) = (" + std::to_string(rec.eta) +
          ", " + std::to_string(rec.phi) + ") in the file, but (" +
          std::to_string(eta) + ", " + std::to_string(phi) +
          ") in the configuration");
    }

    scan_data.tracks.push_back({ray.pos(), 0.f, ray.dir(), 0.f});
  }

  std::clog << "-----------------------------------\n"
            << "Read " << scan_data.tracks.size() << " tracks from "
            << file_name.string() << "\n"
            << "-----------------------------------\n"
            << std::endl;

  return scan_data;
}

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
inline auto material_scan(
    const typename detector_t::geometry_context &gctx, const detector_t &det,
    const typename detector_t::name_map &names,
    const test::material_scan_config<typename detector_t::algebra_type> &cfg)
    -> material_scan_data<typename detector_t::algebra_type> {
  using algebra_t = typename detector_t::algebra_type;
  using point2_t = dpoint2D<algebra_t>;
  using ray_t = detail::ray<algebra_t>;
  using track_generator_t = uniform_track_generator<ray_t>;

  const std::string det_name{det.name(names)};

  if (cfg.truth_from_file()) {
    return read_material_scan(det_name, cfg);
  }

  std::size_t n_tracks{0u};
  auto ray_generator = track_generator_t(cfg.track_generator());

  DETRAY_INFO_HOST("Running material scan on: " << det_name << "\n("
                                                << ray_generator.size()
                                                << " rays) ...\n");

  // Trace material per ray
  material_scan_data<algebra_t> scan_data{};
  scan_data.tracks.reserve(ray_generator.size());
  scan_data.material.reserve(ray_generator.size());
  scan_data.n_steps.reserve(ray_generator.size());

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

    // New material record
    track_material track_mat{};
    track_mat.eta = vector::eta(ray.dir());
    track_mat.phi = vector::phi(ray.dir());
    std::size_t n_steps{0u};

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

      ++n_steps;

      if (mx0 > 0.f) {
        track_mat.sX0 += seg / mx0;
        track_mat.tX0 += t / mx0;
      } else {
        DETRAY_ERROR_HOST("Encountered invalid X_0: " << mx0
                                                      << "\nOn surface: " << sf);
      }
      if (ml0 > 0.f) {
        track_mat.sL0 += seg / ml0;
        track_mat.tL0 += t / ml0;
      } else {
        DETRAY_ERROR_HOST("Encountered invalid L_0: " << ml0
                                                      << "\nOn surface: " << sf);
      }
    }

    if (track_mat.sX0 == 0.f || track_mat.sL0 == 0.f || track_mat.tX0 == 0.f ||
        track_mat.tL0 == 0.f) {
      DETRAY_VERBOSE_HOST("No material recorded for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray);
    }

    scan_data.material.push_back(track_mat);
    scan_data.n_steps.push_back(n_steps);

    ++n_tracks;
  }

  std::clog << "-----------------------------------\n"
            << "Tested " << n_tracks << " tracks\n"
            << "-----------------------------------\n"
            << std::endl;

  // Write recorded material to csv file
  write_material(material_file_path(det_name, cfg.material_file()).string(),
                 scan_data.material);

  return scan_data;
}

/// Trace the material during navigation and compare it to the truth material
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
/// @param scan_data the truth tracks and material
/// @param cfg the material validation configuration
///
/// @returns the number of compared and failed tracks
template <typename backend_t, concepts::detector detector_t>
inline material_validation_result material_validation(
    vecmem::memory_resource *host_mr, const detector_t &det,
    const typename detector_t::name_map &names,
    const material_scan_data<typename detector_t::algebra_type> &scan_data,
    const test::material_validation_config<typename detector_t::algebra_type>
        &cfg) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  const auto &tracks = scan_data.tracks;
  const auto &truth_track_mat = scan_data.material;

  const std::string det_name{det.name(names)};

  DETRAY_INFO_HOST("Running material validation on: " << det_name << "...\n");

  // Capacities of the device material step allocations: Number of material
  // steps found by the material scan plus contingency
  constexpr std::size_t default_capacity{80u};
  std::vector<std::size_t> capacities(tracks.size(), default_capacity);
  if (scan_data.n_steps.size() == tracks.size()) {
    for (std::size_t i = 0u; i < tracks.size(); ++i) {
      capacities[i] = scan_data.n_steps[i] + cfg.step_contingency();
    }
  }

  // Run the propagation and record the accumulated material
  auto [track_mat_vec, mat_steps] = backend_t{}(
      host_mr, cfg.device_mr(), det, cfg.propagation(), tracks, capacities);

  // One material record per track
  if (tracks.size() != track_mat_vec.size()) {
    throw std::runtime_error(
        "Material validation: Number of material records (" +
        std::to_string(track_mat_vec.size()) +
        ") does not match the number of tracks (" +
        std::to_string(tracks.size()) + ")");
  }

  auto get_rel_error = [](const double truth, const double rec) {
    constexpr scalar_t e{std::numeric_limits<double>::epsilon()};

    if (truth <= e && rec <= e) {
      // No material for this ray => valid
      return double{0.f};
    } else if (truth <= e) {
      // Material found where none should be
      return detail::invalid_value<double>();
    } else {
      return math::fabs(truth - rec) / truth;
    }
  };

  // Collect some statistics
  material_validation_result result{};
  const scalar_t rel_error{cfg.relative_error()};
  for (std::size_t i = 0u; i < track_mat_vec.size(); ++i) {
    if (result.n_tracks >= cfg.n_tracks()) {
      break;
    }

    const track_material &truth_mat = truth_track_mat[i];
    const track_material &recorded_mat = track_mat_vec[i];

    bool is_good{true};
    auto check = [&](const double truth, const double rec,
                     const std::string_view what) {
      if (!(get_rel_error(truth, rec) < rel_error)) {
        is_good = false;
        // Same floating point precision as the googletest messages
        std::ostringstream msg{};
        msg << std::setprecision(std::numeric_limits<double>::digits10 + 2)
            << "Track " << result.n_tracks << " (" << what << "): Truth "
            << truth << ", Nav. " << rec;
        std::cerr << msg.str() << std::endl;
      }
    };

    check(truth_mat.sX0, recorded_mat.sX0, "X0 / path");
    check(truth_mat.tX0, recorded_mat.tX0, "X0 / thickness");
    check(truth_mat.sL0, recorded_mat.sL0, "L0 / path");
    check(truth_mat.tL0, recorded_mat.tL0, "L0 / thickness");

    if (!is_good) {
      ++result.n_failed;
    }

    ++result.n_tracks;
  }

  std::clog << "-----------------------------------\n"
            << "Tested " << result.n_tracks << " tracks\n"
            << "-----------------------------------\n"
            << std::endl;

  // Print accumulated material per track
  write_material(
      material_file_path(det_name, cfg.material_file(), backend_t::name)
          .string(),
      track_mat_vec);

  return result;
}

}  // namespace detray::material_validator
