// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/whiteboard.hpp"
#include "detray/test/validation/material_validation_config.hpp"
#include "detray/test/validation/material_validation_utils.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

namespace detray::test {

/// Run the material validation on the host
struct run_material_validation {
  static constexpr std::string_view name{"cpu"};

  template <typename detector_t>
  auto operator()(
      vecmem::memory_resource *host_mr, vecmem::memory_resource * /*mr*/,
      const detector_t &det, const propagation::config &cfg,
      const dvector<free_track_parameters<typename detector_t::algebra_type>>
          &tracks,
      const std::vector<std::size_t> & /*mask*/ = {}) {
    using scalar_t = dscalar<typename detector_t::algebra_type>;

    typename detector_t::geometry_context gctx{};

    dvector<material_validator::track_material<scalar_t>> track_mat_vec{
        host_mr};
    track_mat_vec.reserve(tracks.size());

    dvector<dvector<material_record<scalar_t>>> mat_steps_vec{host_mr};
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

/// @brief Test class that runs the material validation for a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed outside this class
template <typename detector_t, typename material_validator_t>
class material_validation_impl : public test::fixture_base<> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using free_track_parameters_t = free_track_parameters<algebra_t>;
  using track_material_t = material_validator::track_material<scalar_t>;

 public:
  using fixture_type = test::fixture_base<>;
  using config = detray::test::material_validation_config<algebra_t>;

  explicit material_validation_impl(
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

    // Name of the material scan data collection
    m_scan_data_name = m_det.name(m_names) + "_material_scan";
    m_track_data_name = m_det.name(m_names) + "_material_scan_tracks";

    // Check that data is available in memory
    if (!m_whiteboard->exists(m_scan_data_name)) {
      throw std::invalid_argument(
          "Material validation: Could not find scan data on whiteboard."
          "Please run material scan first.");
    }
    if (!m_whiteboard->exists(m_track_data_name)) {
      throw std::invalid_argument(
          "Material validation: Could not find track data on whiteboard."
          "Please run material scan first.");
    }
  }

  /// Run the check
  void TestBody() override {
    using namespace detray;

    // Fetch the input data
    const auto &tracks =
        m_whiteboard->template get<dvector<free_track_parameters_t>>(
            m_track_data_name);

    const auto &truth_track_mat =
        m_whiteboard->template get<dvector<track_material_t>>(m_scan_data_name);

    DETRAY_INFO_HOST("Running material validation on: " << m_det.name(m_names)
                                                        << "...\n");

    // only needed for device material steps allocations
    // @TODO: For now, guess how many surface might be encountered
    std::vector<std::size_t> capacities(tracks.size(), 80u);

    // Run the propagation on device and record the accumulated material
    auto [track_mat_vec, mat_steps] =
        material_validator_t{}(&m_host_mr, m_cfg.device_mr(), m_det,
                               m_cfg.propagation(), tracks, capacities);

    // One material record per track
    ASSERT_EQ(tracks.size(), track_mat_vec.size());

    // Collect some statistics
    std::size_t n_tracks{0u};
    const scalar_t rel_error{m_cfg.relative_error()};
    for (std::size_t i = 0u; i < track_mat_vec.size(); ++i) {
      if (n_tracks >= m_cfg.n_tracks()) {
        break;
      }

      const track_material_t &truth_mat = truth_track_mat[i];
      const track_material_t &recorded_mat = track_mat_vec[i];

      auto get_rel_error = [](const scalar_t truth, const scalar_t rec) {
        constexpr scalar_t e{std::numeric_limits<scalar_t>::epsilon()};

        if (truth <= e && rec <= e) {
          // No material for this ray => valid
          return scalar_t{0.f};
        } else if (truth <= e) {
          // Material found where none should be
          return detail::invalid_value<scalar_t>();
        } else {
          return math::fabs(truth - rec) / truth;
        }
      };

      EXPECT_LT(get_rel_error(truth_mat.sX0, recorded_mat.sX0), rel_error)
          << "Track " << n_tracks << " (X0 / path): Truth " << truth_mat.sX0
          << ", Nav. " << recorded_mat.sX0;
      EXPECT_LT(get_rel_error(truth_mat.tX0, recorded_mat.tX0), rel_error)
          << "Track " << n_tracks << " (X0 / thickness): Truth "
          << truth_mat.tX0 << ", Nav. " << recorded_mat.tX0;
      EXPECT_LT(get_rel_error(truth_mat.sL0, recorded_mat.sL0), rel_error)
          << "Track " << n_tracks << " (L0 / path): Truth " << truth_mat.sL0
          << ", Nav. " << recorded_mat.sL0;
      EXPECT_LT(get_rel_error(truth_mat.tL0, recorded_mat.tL0), rel_error)
          << "Track " << n_tracks << " (L0 / thickness): Truth "
          << truth_mat.tL0 << ", Nav. " << recorded_mat.tL0;

      ++n_tracks;
    }

    std::clog << "-----------------------------------\n"
              << "Tested " << n_tracks << " tracks\n"
              << "-----------------------------------\n"
              << std::endl;

    // Print accumulated material per track
    std::filesystem::path mat_path{m_cfg.material_file()};
    const auto data_path{mat_path.parent_path()};

    auto file_name{data_path /
                   (m_det.name(m_names) + "_" + mat_path.stem().string() + "_" +
                    std::string(material_validator_t::name) + ".csv")};

    material_validator::write_material(file_name.string(), track_mat_vec);
  }

 private:
  /// The configuration of this test
  config m_cfg;
  /// Vecmem memory resource for the host allocations
  vecmem::host_memory_resource m_host_mr{};
  /// Name of the input data collections
  std::string m_scan_data_name{""};
  std::string m_track_data_name{""};
  /// The geometry context to check
  typename detector_t::geometry_context m_gctx{};
  /// The detector to be checked
  const detector_t &m_det;
  /// Volume names
  const typename detector_t::name_map &m_names;
  /// Whiteboard to pin data
  std::shared_ptr<test::whiteboard> m_whiteboard{nullptr};
};

template <typename detector_t>
using material_validation =
    material_validation_impl<detector_t, run_material_validation>;

}  // namespace detray::test
