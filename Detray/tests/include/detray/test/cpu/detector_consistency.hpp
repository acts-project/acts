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
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/utils/hash_tree.hpp"

// System include(s)
#include <iostream>
#include <string>
#include <string_view>

namespace detray::test {

struct consistency_check_config
    : public detray::test::fixture_base<>::configuration {
  std::string m_name{"detector_consistency"};
  bool m_write_graph{false};

  /// Getters
  /// @{
  const std::string &name() const { return m_name; }
  bool write_graph() const { return m_write_graph; }
  /// @}

  /// Setters
  /// @{
  consistency_check_config &name(const std::string_view v) {
    m_name = v;
    return *this;
  }
  consistency_check_config &write_graph(const bool do_write) {
    m_write_graph = do_write;
    return *this;
  }
  /// @}
};

/// @brief Test class that runs the consistency check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class consistency_check : public detray::test::fixture_base<> {
 public:
  using fixture_type = detray::test::fixture_base<>;
  using config = consistency_check_config;

  template <typename config_t>
  explicit consistency_check(
      const detector_t &det, const typename detector_t::name_map &names,
      const config_t &cfg = {},
      const typename detector_t::geometry_context &gctx = {})
      : m_cfg{cfg}, m_gctx{gctx}, m_det{det}, m_names{names} {}

  /// Run the consistency check
  void TestBody() override {
    DETRAY_INFO_HOST("Running consistency check on: " << m_det.name(m_names)
                                                      << std::endl);

    // Build the graph
    volume_graph graph(m_det);

    ASSERT_TRUE(detail::check_consistency(m_det, true, m_names))
        << graph.to_string();

    if (m_cfg.write_graph()) {
      DETRAY_INFO_HOST(graph.to_string());
    }

    // Not currently supported
    if (false) {
      constexpr std::size_t root_hash{3244u};  // < toy detector only

      const auto &adj_mat = graph.adjacency_matrix();

      auto geo_checker = hash_tree(adj_mat);

      EXPECT_EQ(geo_checker.root(), root_hash) << graph.to_string();
    }
  }

 private:
  /// The configuration of this test
  const config m_cfg;
  /// The geometry context to check
  typename detector_t::geometry_context m_gctx{};
  /// The detector to be checked
  const detector_t &m_det;
  /// Volume names
  const typename detector_t::name_map &m_names;
};

}  // namespace detray::test
