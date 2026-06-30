// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/whiteboard.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/material_validation_utils.hpp"
#include "detray/test/validation/navigation_validation_config.hpp"
#include "detray/test/validation/navigation_validation_utils.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <memory>
#include <optional>
#include <string>

namespace detray::test {

/// @brief Test class that runs the navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t, template <typename> class scan_type>
class navigation_validation : public test::fixture_base<> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;
  using free_track_parameters_t = free_track_parameters<algebra_t>;
  using trajectory_type = typename scan_type<algebra_t>::trajectory_type;
  using truth_trace_t = typename scan_type<
      algebra_t>::template intersection_trace_type<detector_t>;

  /// Switch between rays and helices
  static constexpr auto k_use_rays{
      std::is_same_v<detail::ray<algebra_t>, trajectory_type>};

 public:
  using fixture_type = test::fixture_base<>;
  using config = navigation_validation_config<algebra_t>;

  explicit navigation_validation(
      const detector_t &det, const typename detector_t::name_map &names,
      const config &cfg = {}, std::shared_ptr<test::whiteboard> wb = nullptr,
      const typename detector_t::geometry_context gctx = {})
      : m_cfg{cfg},
        m_gctx{gctx},
        m_det{det},
        m_names{names},
        m_whiteboard{std::move(wb)} {
    if (!m_whiteboard) {
      throw std::invalid_argument("No white board was passed to " +
                                  m_cfg.name() + " test");
    }
  }

  /// Run the check
  void TestBody() override {
    using namespace detray;
    using namespace navigation;

    using intersection_t =
        typename truth_trace_t::value_type::intersection_type;

    // Runge-Kutta stepper
    using hom_bfield_t = typename bfield::const_field_t<scalar_t>::view_t;
    using bfield_t =
        std::conditional_t<k_use_rays, navigation_validator::empty_bfield,
                           hom_bfield_t>;
    using rk_stepper_t =
        rk_stepper<hom_bfield_t, algebra_t, unconstrained_step<scalar_t>,
                   stepper_rk_policy<scalar_t>, stepping::print_inspector>;
    using line_stepper_t = line_stepper<algebra_t, unconstrained_step<scalar_t>,
                                        stepper_default_policy<scalar_t>,
                                        stepping::print_inspector>;
    using stepper_t =
        std::conditional_t<k_use_rays, line_stepper_t, rk_stepper_t>;

    std::optional<bfield_t> b_field{};
    if constexpr (!k_use_rays) {
      b_field.emplace(create_const_field<scalar_t>(m_cfg.B_vector()));
    }

    // Use ray or helix
    const std::string det_name{m_det.name(m_names)};
    const std::string truth_data_name{k_use_rays ? det_name + "_ray_scan"
                                                 : det_name + "_helix_scan"};

    // Collect some statistics
    std::size_t n_tracks{0u};
    std::size_t n_matching_error{0u};
    std::size_t n_fatal{0u};
    // Total number of encountered surfaces
    navigation_validator::surface_stats n_surfaces{};
    // Missed by navigator
    navigation_validator::surface_stats n_miss_nav{};
    // Missed by truth finder
    navigation_validator::surface_stats n_miss_truth{};

    DETRAY_INFO_HOST("Fetching data from white board...");
    if (!m_whiteboard->exists(truth_data_name)) {
      throw std::runtime_error(
          "White board is empty! Please run detector scan first");
    }
    auto &truth_traces =
        m_whiteboard->template get<std::vector<truth_trace_t>>(truth_data_name);
    ASSERT_EQ(m_cfg.n_tracks(), truth_traces.size());

    DETRAY_INFO_HOST("Running navigation validation on: " << det_name
                                                          << "...\n");

    std::string momentum_str{""};
    const std::string prefix{k_use_rays ? det_name + "_ray_"
                                        : det_name + "_helix_"};

    const auto data_path{
        std::filesystem::path{m_cfg.track_param_file()}.parent_path()};

    // Create an output file path
    auto make_path = [&data_path, &prefix, &momentum_str](
                         const std::string &name,
                         const std::string &extension = ".csv") {
      return data_path / (prefix + name + momentum_str + extension);
    };

    std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
    const std::string debug_file_name{
        make_path("navigation_validation", ".txt")};
    detray::io::file_handle debug_file{debug_file_name, io_mode};

    // Keep a record of track positions and material along the track
    dvector<dvector<intersection_record<detector_t>>> recorded_traces{};
    dvector<material_validator::track_material<scalar_t>> track_mat_vec{};
    std::vector<std::pair<trajectory_type, std::vector<intersection_t>>>
        missed_intersections{};

    scalar_t min_pT{std::numeric_limits<scalar_t>::max()};
    scalar_t max_pT{-std::numeric_limits<scalar_t>::max()};
    for (auto &truth_trace : truth_traces) {
      if (n_tracks >= m_cfg.n_tracks()) {
        break;
      }

      // Follow the test trajectory with a track and check, if we find
      // the same volumes and distances along the way
      const auto &start = truth_trace.front();
      const auto &track = start.track_param();
      assert(!track.is_invalid());
      trajectory_type test_traj = get_parametrized_trajectory(track);

      const scalar q = start.charge;
      const scalar pT{q == 0.f ? 1.f * unit<scalar>::GeV : track.pT(q)};
      const scalar p{q == 0.f ? 1.f * unit<scalar>::GeV : track.p(q)};

      // If the momentum is unknown, 1 GeV is the safest option to keep
      // the direction vector normalized
      if (detray::detail::is_invalid_value(m_cfg.p_range()[0])) {
        min_pT = std::min(min_pT, pT);
        max_pT = std::max(max_pT, pT);
      } else {
        min_pT = m_cfg.p_range()[0];
        max_pT = m_cfg.p_range()[1];
      }
      assert(min_pT > 0.f);
      assert(max_pT > 0.f);
      assert(min_pT < std::numeric_limits<scalar_t>::max());
      assert(max_pT < std::numeric_limits<scalar_t>::max());

      // Run the propagation
      auto [success, obj_tracer, step_trace, mat_record, mat_trace, nav_printer,
            step_printer] =
          navigation_validator::record_propagation<stepper_t>(
              m_gctx, &m_host_mr, m_det, m_cfg.propagation(), track,
              m_cfg.ptc_hypothesis(), b_field);

      if (success) {
        assert(!obj_tracer.object_trace.empty());
        // The navigator does not record the initial track position:
        // add it as a dummy record
        obj_tracer.object_trace.insert(
            obj_tracer.object_trace.begin(),
            {track.pos(), track.dir(), start.intersection});

        // Adjust the track charge, which is unknown to the navigation
        for (auto &record : obj_tracer.object_trace) {
          record.charge = q;
          record.p_mag = p;
        }

        auto [result, n_missed_nav, n_missed_truth, n_error, missed_inters] =
            navigation_validator::compare_traces(
                m_cfg, truth_trace, obj_tracer.object_trace, test_traj,
                n_tracks, &(*debug_file));

        missed_intersections.push_back(
            std::make_pair(test_traj, std::move(missed_inters)));

        // Update statistics
        success = success && result;
        n_miss_nav += n_missed_nav;
        n_miss_truth += n_missed_truth;
        n_matching_error += n_error;

      } else {
        // Propagation did not succeed
        ++n_fatal;

        std::vector<intersection_t> missed_inters{};
        missed_intersections.push_back(
            std::make_pair(test_traj, missed_inters));
      }

      if (!success) {
        // Write debug info to file
        *debug_file << "TEST TRACK " << n_tracks << ":\n\n"
                    << "NAVIGATOR\n\n"
                    << nav_printer.to_string() << "\nSTEPPER\n\n"
                    << step_printer.to_string();

        detector_scanner::display_error(
            m_gctx, m_det, m_names, m_cfg.name(), test_traj, truth_trace,
            m_cfg.svg_style(), n_tracks, m_cfg.n_tracks(),
            obj_tracer.object_trace);
      }

      recorded_traces.push_back(std::move(obj_tracer.object_trace));
      track_mat_vec.push_back(mat_record);

      EXPECT_TRUE(success)
          << "\nDETRAY INFO (HOST): Wrote navigation debugging data in: "
          << debug_file_name;

      ++n_tracks;

      // After dummy records insertion, traces should have the same size
      ASSERT_EQ(truth_trace.size(), recorded_traces.back().size());

      // Count the number of different surface types on this trace
      navigation_validator::surface_stats n_truth{};
      navigation_validator::surface_stats n_nav{};
      for (std::size_t i = 0; i < truth_trace.size(); ++i) {
        const auto truth_desc = truth_trace[i].intersection.surface();
        const auto rec_desc = recorded_traces.back()[i].intersection.surface();

        // Exclude dummy records for missing surfaces
        if (!truth_desc.identifier().is_invalid()) {
          n_truth.count(truth_desc);
        }
        if (!rec_desc.identifier().is_invalid()) {
          n_nav.count(rec_desc);
        }
      }

      // Take max count, since either trace might have skipped surfaces
      const std::size_t n_portals{
          math::max(n_truth.n_portals, n_nav.n_portals)};
      const std::size_t n_sensitives{
          math::max(n_truth.n_sensitives, n_nav.n_sensitives)};
      const std::size_t n_passives{
          math::max(n_truth.n_passives, n_nav.n_passives)};
      const std::size_t n{n_portals + n_sensitives + n_passives};

      // Cannot have less surfaces than truth intersections after matching
      // (Don't count first entry, which records the initial track params)
      ASSERT_TRUE(n >= (truth_trace.size() - 1u));

      n_surfaces.n_portals += n_portals;
      n_surfaces.n_sensitives += n_sensitives;
      n_surfaces.n_passives += n_passives;
    }

    // Calculate and display the result
    navigation_validator::print_efficiency(n_tracks, n_surfaces, n_miss_nav,
                                           n_miss_truth, n_fatal,
                                           n_matching_error);

    // Print track positions for plotting
    if constexpr (!k_use_rays) {
      momentum_str =
          "_" +
          std::to_string(std::floor(10. * static_cast<double>(min_pT)) / 10.) +
          "_" +
          std::to_string(std::ceil(10. * static_cast<double>(max_pT)) / 10.) +
          "_GeV";
    }

    const auto truth_trk_path{make_path("truth_track_params")};
    const auto trk_path{make_path("navigation_track_params")};
    const auto truth_intr_path{make_path("truth_intersections")};
    const auto intr_path{make_path("navigation_intersections")};
    const auto mat_path{make_path("accumulated_material")};
    const auto missed_path{make_path("missed_intersections_dists")};

    // Write the distance of the missed intersection local position
    // to the surface boundaries to file for plotting
    navigation_validator::write_dist_to_boundary(
        m_det, m_names, missed_path.string(), missed_intersections);
    detector_scanner::write_tracks(truth_trk_path.string(), truth_traces);
    navigation_validator::write_tracks(trk_path.string(), recorded_traces);
    detector_scanner::write_intersections(truth_intr_path.string(),
                                          truth_traces);
    detector_scanner::write_intersections(intr_path.string(), recorded_traces);
    material_validator::write_material(mat_path.string(), track_mat_vec);

    DETRAY_INFO_HOST("Wrote distance to boundary of missed intersections in: "
                     << missed_path);
    DETRAY_INFO_HOST("Wrote truth track states in: " << truth_trk_path);
    DETRAY_INFO_HOST("Wrote recorded track states in: " << trk_path);
    DETRAY_INFO_HOST(
        "Wrote recorded truth intersections in: " << truth_intr_path);
    DETRAY_INFO_HOST("Wrote recorded track intersections in: " << intr_path);
    DETRAY_INFO_HOST("Wrote accumulated material in: " << mat_path);
  }

 private:
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
  /// The geometry context to check
  typename detector_t::geometry_context m_gctx{};
  /// Vecmem memory resource for the host allocations
  vecmem::host_memory_resource m_host_mr{};
  /// The detector to be checked
  const detector_t &m_det;
  /// Volume names
  const typename detector_t::name_map &m_names;
  /// Whiteboard to pin data
  std::shared_ptr<test::whiteboard> m_whiteboard{nullptr};
};

template <typename detector_t>
using straight_line_navigation =
    navigation_validation<detector_t, detray::ray_scan>;

template <typename detector_t>
using helix_navigation = navigation_validation<detector_t, detray::helix_scan>;

}  // namespace detray::test
