// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/navigation/volume_graph.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"

// Detray IO include(s)
#include "detray/io/utils/create_path.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/framework/whiteboard.hpp"
#include "detray/test/validation/detector_scan_config.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"

// System include(s)
#include <iostream>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Test class that runs the ray/helix scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t, template <typename> class scan_type>
class detector_scan : public test::fixture_base<> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using free_track_parameters_t = free_track_parameters<algebra_t>;
  using intersection_trace_t = typename scan_type<
      algebra_t>::template intersection_trace_type<detector_t>;
  using trajectory_type = typename scan_type<algebra_t>::trajectory_type;
  using uniform_gen_t =
      detail::random_numbers<scalar_t,
                             std::uniform_real_distribution<scalar_t>>;
  using track_generator_t =
      random_track_generator<free_track_parameters_t, uniform_gen_t>;

  /// Switch between rays and helices
  static constexpr auto k_use_rays{
      std::is_same_v<detail::ray<algebra_t>, trajectory_type>};

 public:
  using fixture_type = test::fixture_base<>;
  using config = detector_scan_config<track_generator_t, algebra_t>;

  explicit detector_scan(const detector_t &det,
                         const typename detector_t::name_map &names,
                         const config &cfg = {},
                         std::shared_ptr<test::whiteboard> wb = nullptr,
                         const typename detector_t::geometry_context gctx = {})
      : m_cfg{cfg},
        m_gctx{gctx},
        m_det{det},
        m_names{names},
        m_whiteboard{std::move(wb)} {
    if (!m_whiteboard) {
      throw std::invalid_argument("No white board was passed to " +
                                  m_cfg.name());
    }
  }

  /// Run the detector scan
  void TestBody() override {
    // Get the volume adjaceny matrix from the detector
    volume_graph graph(m_det);
    const auto &adj_mat = graph.adjacency_matrix();

    // Fill adjacency matrix from ray scan and compare
    dvector<dindex> adj_mat_scan(adj_mat.size(), 0);
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    // Index of the volume that the test trajectory origin lies in
    dindex start_index{0u};

    DETRAY_INFO_HOST("Running scan on: " << m_det.name(m_names) << "\n");

    // Fill detector scan data to white board
    const std::size_t n_helices = fill_scan_data();

    auto &detector_scan_traces =
        m_whiteboard->template get<std::vector<intersection_trace_t>>(
            m_cfg.name());

    const std::string det_name{m_det.name(m_names)};
    const std::string prefix{k_use_rays ? det_name + "_ray"
                                        : det_name + "_helix"};

    const auto data_path{
        std::filesystem::path{m_cfg.track_param_file()}.parent_path()};

    // Make sure the output directories exit
    io::create_path(data_path);

    std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
    detray::io::file_handle debug_file{
        data_path / (prefix + "_detector_scan.txt"), io_mode};

    DETRAY_INFO_HOST("Checking trace data...\n");

    // Iterate through the scan data and perform checks
    std::size_t n_tracks{0u};
    for (int i = static_cast<int>(detector_scan_traces.size()) - 1; i >= 0;
         --i) {
      const auto j{static_cast<std::size_t>(i)};

      auto &intersection_trace = detector_scan_traces[j];
      assert((intersection_trace.size() > 0) && "Invalid intersection trace");

      // Retrieve the test trajectory
      const auto &trck_param = intersection_trace.front().track_param();
      trajectory_type test_traj = get_parametrized_trajectory(trck_param);

      // Run overlaps check on the trace and remove certain
      // allowed duplication (oversized portals from ACTS)
      if (m_cfg.overlaps_removal()) {
        const dindex_range overlap_idx = detector_scanner::overlaps_removal(
            intersection_trace, m_cfg.overlaps_tol());

        // Drop an svg of the trajectory where the overlap was found
        if (overlap_idx[1] - overlap_idx[0] != 0u) {
          constexpr bool verbose{false};
          detector_scanner::display_error(
              m_gctx, m_det, m_names, "OVERLAP_" + m_cfg.name(), test_traj,
              intersection_trace, m_cfg.svg_style(), n_tracks, n_helices,
              intersection_trace_t{}, overlap_idx, verbose);
        }
      }

      // Run consistency checks on the trace
      bool success = detector_scanner::check_trace<detector_t>(
          intersection_trace, start_index, adj_mat_scan, obj_hashes);

      // Display the detector, track and intersections for debugging
      if (!success) {
        detector_scanner::display_error(
            m_gctx, m_det, m_names, m_cfg.name(), test_traj, intersection_trace,
            m_cfg.svg_style(), n_tracks, n_helices, intersection_trace_t{});
      }

      EXPECT_TRUE(success);

      // Remove faulty trace for the following steps
      if (!success) {
        *debug_file << detector_scanner::print_trace(intersection_trace, j);

        DETRAY_ERROR_HOST("Skipped faulty trace no. " << j);
        detector_scan_traces.erase(detector_scan_traces.begin() + i);
      } else {
        ++n_tracks;
      }
    }
    std::clog << "------------------------------------\n"
              << "Tested " << n_tracks << " tracks: OK\n"
              << "------------------------------------\n"
              << std::endl;

    // Check that the links that were discovered by the scan match the
    // volume graph
    // ASSERT_TRUE(adj_mat == adj_mat_scan) <<
    // detector_scanner::print_adj(adj_mat_scan);

    // Compare the adjacency that was discovered in the ray scan to the
    // hashed one for the toy detector.
    // The hash tree is still Work in Progress !
    /*auto geo_checker = hash_tree(adj_mat);
    const bool check_links = (geo_checker.root() == root_hash);

    std::clog << "All links reachable: " << (check_links ? "OK" : "FAILURE")
            << std::endl;*/
  }

 private:
  /// Get the truth data, either from file or by generating it
  /// @returns the number of helices
  std::size_t fill_scan_data() {
    /// Record the scan for later use
    std::vector<intersection_trace_t> intersection_traces;

    auto trk_state_generator = track_generator_t(m_cfg.track_generator());
    const std::size_t n_helices{trk_state_generator.size()};
    intersection_traces.reserve(n_helices);

    std::string momentum_str{""};
    if constexpr (!k_use_rays) {
      const auto pT_range = m_cfg.track_generator().mom_range();
      // Remove floating point imprecisions
      momentum_str =
          "_" +
          std::to_string(std::floor(10. * static_cast<double>(pT_range[0])) /
                         10.) +
          "_" +
          std::to_string(std::ceil(10. * static_cast<double>(pT_range[1])) /
                         10.) +
          "_GeV";
    }

    std::string track_param_file_name{m_cfg.track_param_file() + momentum_str +
                                      ".csv"};

    std::string intersection_file_name{m_cfg.intersection_file() +
                                       momentum_str + ".csv"};

    const bool data_files_exist{io::file_exists(intersection_file_name) &&
                                io::file_exists(track_param_file_name)};

    if (data_files_exist) {
      DETRAY_INFO_HOST("Reading data from file...");

      // Fill the intersection traces from file
      detector_scanner::read(intersection_file_name, track_param_file_name,
                             intersection_traces);
    } else {
      DETRAY_INFO_HOST("Generating trace data...");

      for (auto trk : trk_state_generator) {
        // Get ground truth from track
        trajectory_type test_traj = get_parametrized_trajectory(trk);

        // The track generator can randomize the sign of the charge
        const scalar_t qabs{math::fabs(m_cfg.m_trk_gen_cfg.charge())};
        const scalar_t q{math::copysign(qabs, trk.qop())};

        // Shoot trajectory through the detector and record all
        // surfaces it encounters
        // @note: For rays, set the momentum to 1 GeV to keep the
        //        direction vector normalized
        const scalar p{q == 0.f ? 1.f * unit<scalar>::GeV : trk.p(q)};
        auto trace = detector_scanner::run<scan_type>(
            m_gctx, m_det, test_traj, m_cfg.mask_tolerance(), p);

        intersection_traces.push_back(std::move(trace));
      }
    }

    // Save the results

    // Csv output
    if (!data_files_exist && m_cfg.write_intersections()) {
      detector_scanner::write_tracks(track_param_file_name,
                                     intersection_traces);
      detector_scanner::write_intersections(intersection_file_name,
                                            intersection_traces);

      DETRAY_INFO_HOST("  ->Wrote  " << intersection_traces.size()
                                     << " intersection traces to file");
    }

    DETRAY_INFO_HOST("  ->Adding " << intersection_traces.size()
                                   << " intersection traces to whiteboard");

    // Move the data to the whiteboard
    m_whiteboard->add(m_cfg.name(), std::move(intersection_traces));

    return n_helices;
  }

  /// @returns either the helix or ray corresponding to the input track
  /// parameters @param track
  trajectory_type get_parametrized_trajectory(
      const free_track_parameters_t &track) {
    std::unique_ptr<trajectory_type> test_traj{nullptr};
    if constexpr (k_use_rays) {
      test_traj = std::make_unique<trajectory_type>(track);
    } else {
      test_traj = std::make_unique<trajectory_type>(track, m_cfg.B_vector());
    }
    return *(test_traj.release());
  }

  /// The configuration of this test
  config m_cfg;
  /// The geometry context to scan
  typename detector_t::geometry_context m_gctx{};
  /// The detector to be checked
  const detector_t &m_det;
  /// Volume names
  const typename detector_t::name_map &m_names;
  /// Whiteboard to pin data
  std::shared_ptr<test::whiteboard> m_whiteboard{nullptr};
};

template <typename detector_t>
using ray_scan = detector_scan<detector_t, detray::ray_scan>;

template <typename detector_t>
using helix_scan = detector_scan<detector_t, detray::helix_scan>;

}  // namespace detray::test
