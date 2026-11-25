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
#include "detray/test/utils/data_record.hpp"
#include "detray/test/validation/material_scanner.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <string>

namespace detray::test {

/// @brief Test class that runs the material ray scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <concepts::detector detector_t>
class material_scan_fixture : public test::fixture_base<> {
  using algebra_t = typename detector_t::algebra_type;

 public:
  /// Test fixture base type
  using fixture_type = test::fixture_base<>;
  /// Material scan configuration type
  using config = validation::material_scan_config<algebra_t>;

  /// Parameterized constructor
  /// @param det the detector to be scanned
  /// @param name name map of the detector
  /// @param cfg the scan config
  /// @param wb pointer to the whiteboard the data will be pinned to
  /// @param gctx geometry context
  explicit material_scan_fixture(
      const detector_t &det, const typename detector_t::name_map &names,
      const config &cfg = {}, std::shared_ptr<test::whiteboard> wb = nullptr,
      const typename detector_t::geometry_context &gctx = {})
      : m_cfg{cfg},
        m_gctx{gctx},
        m_det{det},
        m_names{names},
        m_whiteboard{std::move(wb)} {
    if (!m_whiteboard) {
      throw std::invalid_argument("No white board was passed to " +
                                  m_cfg.name() + "(" + m_det.name(m_names) +
                                  ")");
    }
  }

  /// Run the ray scan
  void TestBody() override {
    material_trace_data<algebra_t> scan_data =
        validation::material_scanner::run(m_gctx, m_det, m_names, m_cfg);

    // Pin data to whiteboard
    m_whiteboard->add(m_det.name(m_names) + "_material_scan",
                      std::move(scan_data));
  }

 private:
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

}  // namespace detray::test
