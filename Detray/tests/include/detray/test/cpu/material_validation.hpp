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

// Detray test include(s)
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/whiteboard.hpp"
#include "detray/test/validation/material_validation.hpp"
#include "detray/test/validation/material_validation_config.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <memory>
#include <stdexcept>
#include <string>

namespace detray::test {

/// @brief Test class that runs the material validation for a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed outside this class
template <concepts::detector detector_t, typename material_validator_t>
class material_validation_impl : public test::fixture_base<> {
  using algebra_t = typename detector_t::algebra_type;

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

    // Check that data is available in memory
    if (!m_whiteboard->exists(m_scan_data_name)) {
      throw std::invalid_argument(
          "Material validation: Could not find scan data on whiteboard."
          "Please run material scan first.");
    }
  }

  /// Run the check
  void TestBody() override {
    // Fetch the input data
    const auto &scan_data = m_whiteboard->template get<
        material_validator::material_scan_data<algebra_t>>(m_scan_data_name);

    const auto result =
        material_validator::material_validation<material_validator_t>(
            &m_host_mr, m_det, m_names, scan_data, m_cfg);

    EXPECT_TRUE(result.passed())
        << result.n_failed << " of " << result.n_tracks
        << " tracks have a material discrepancy";
  }

 private:
  /// The configuration of this test
  config m_cfg;
  /// Vecmem memory resource for the host allocations
  vecmem::host_memory_resource m_host_mr{};
  /// Name of the input data collection
  std::string m_scan_data_name{""};
  /// The geometry context to check
  typename detector_t::geometry_context m_gctx{};
  /// The detector to be checked
  const detector_t &m_det;
  /// Volume names
  const typename detector_t::name_map &m_names;
  /// Whiteboard to pin data
  std::shared_ptr<test::whiteboard> m_whiteboard{nullptr};
};

template <concepts::detector detector_t>
using material_validation =
    material_validation_impl<detector_t, run_material_validation>;

}  // namespace detray::test
