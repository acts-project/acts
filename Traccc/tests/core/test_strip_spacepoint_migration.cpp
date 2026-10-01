// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <limits>

#include <gtest/gtest.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

#include "traccc/seeding/device/count_strip_pairs.hpp"
#include "traccc/seeding/device/find_strip_pairs.hpp"
#include "traccc/seeding/device/form_strip_spacepoints_from_pairs.hpp"

namespace {
using namespace traccc;
TEST(strip_spacepoint_migration, count_fill_and_form_standard_and_overlap) {
  vecmem::host_memory_resource mr;
  vecmem::copy copy;
  edm::measurement_collection::host measurements{mr};
  strip_measurement_surface_info_collection_types::host surfaces{&mr};
  strip_pairing_rule_collection_types::host rules{&mr};
  for (unsigned int i = 0; i < 3; ++i) {
    detray::geometry::identifier link{0u};
    link.set_index(i);
    measurements.push_back(
        {{0.f, 0.f}, {0.f, 0.f}, 1u, 0.f, 0.f, 0u, link, {1u, 0u}, i});
    strip_measurement_surface_info info{};
    info.surface_link = link.value();
    info.origin = {100.f, 0.f, 0.f};
    info.local_u = i == 0 ? std::array<scalar, 3>{0.f, 0.f, 1.f}
                          : std::array<scalar, 3>{0.f, 1.f, 0.f};
    info.local_v = i == 0 ? std::array<scalar, 3>{0.f, 1.f, 0.f}
                          : std::array<scalar, 3>{0.f, 0.f, 1.f};
    info.linear = {1u, 1u, 1.f, 20.f, 0.f};
    info.index_mapping = {0u, 1u, 1.f, 1.f, 0.f, true};
    info.spacepoint_variance_r = 0.1f;
    info.spacepoint_variance_z = 0.2f;
    surfaces.push_back(info);
    if (i != 0) {
      strip_pairing_rule rule{};
      rule.reference_surface_link = surfaces[0].surface_link;
      rule.candidate_surface_link = link.value();
      rule.category =
          i == 1 ? strip_pair_category::standard : strip_pair_category::overlap;
      rules.push_back(rule);
    }
  }
  const auto meas_view = vecmem::get_data(measurements);
  const auto surfaces_view = vecmem::get_data(surfaces);
  const auto rules_view = vecmem::get_data(rules);
  unsigned int n_standard = 0, n_overlap = 0;
  // Include an out-of-range lane as in a rounded-up GPU launch.
  for (unsigned int i = 0; i < 4; ++i) {
    device::count_strip_pairs(i, meas_view, surfaces_view, rules_view, point3{},
                              n_standard, n_overlap);
  }
  ASSERT_EQ(n_standard, 1u);
  ASSERT_EQ(n_overlap, 1u);
  strip_pair_collection_types::host standard(n_standard, &mr);
  strip_pair_collection_types::host overlap(n_overlap, &mr);
  unsigned int pos_standard = 0, pos_overlap = 0;
  for (unsigned int i = 0; i < 4; ++i) {
    device::find_strip_pairs(
        i, meas_view, surfaces_view, rules_view, point3{}, pos_standard,
        pos_overlap, vecmem::get_data(standard), vecmem::get_data(overlap));
  }
  ASSERT_EQ(pos_standard, n_standard);
  ASSERT_EQ(pos_overlap, n_overlap);
  EXPECT_EQ(standard[0].measurement_index_2, 1u);
  EXPECT_EQ(overlap[0].measurement_index_2, 2u);
  for (const auto* pairs : {&standard, &overlap}) {
    strip_pair_collection_types::host pairs_with_sentinel{&mr};
    pairs_with_sentinel.push_back(pairs->at(0));
    pairs_with_sentinel.push_back({std::numeric_limits<unsigned int>::max(),
                                   std::numeric_limits<unsigned int>::max(),
                                   0.f, 0.f});
    edm::spacepoint_collection::buffer output(
        2u, mr, vecmem::data::buffer_type::resizable);
    copy.setup(output)->wait();
    // Process a valid pair, an invalid sentinel, and an out-of-range lane.
    for (unsigned int i = 0; i < 3; ++i) {
      device::form_strip_spacepoints_from_pairs(
          i, meas_view, vecmem::get_data(pairs_with_sentinel), surfaces_view,
          point3{}, output);
    }
    ASSERT_EQ(copy.get_size(output), 1u);
    edm::spacepoint_collection::device spacepoints{output};
    const auto sp = spacepoints.at(0);
    EXPECT_FLOAT_EQ(sp.x(), 100.f);
    EXPECT_FLOAT_EQ(sp.y(), 0.f);
    EXPECT_FLOAT_EQ(sp.z(), 0.f);
    EXPECT_FLOAT_EQ(sp.radius_variance(), 0.1f);
    EXPECT_FLOAT_EQ(sp.z_variance(), 0.2f);
    EXPECT_EQ(sp.measurement_index_1(), 0u);
    EXPECT_EQ(sp.measurement_index_2(), pairs->at(0).measurement_index_2);
  }
}
}  // namespace
