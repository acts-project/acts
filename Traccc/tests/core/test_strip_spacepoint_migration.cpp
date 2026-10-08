// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <limits>
#include <numeric>

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
  vecmem::vector<unsigned int> standard_offsets(measurements.size(), 0u, &mr);
  vecmem::vector<unsigned int> overlap_offsets(measurements.size(), 0u, &mr);
  // Include an out-of-range lane as in a rounded-up GPU launch.
  for (unsigned int i = 0; i < 4; ++i) {
    device::count_strip_pairs(i, meas_view, surfaces_view, rules_view, point3{},
                              vecmem::get_data(standard_offsets),
                              vecmem::get_data(overlap_offsets));
  }
  std::partial_sum(standard_offsets.begin(), standard_offsets.end(),
                   standard_offsets.begin());
  std::partial_sum(overlap_offsets.begin(), overlap_offsets.end(),
                   overlap_offsets.begin());
  const unsigned int n_standard = standard_offsets.back();
  const unsigned int n_overlap = overlap_offsets.back();
  ASSERT_EQ(n_standard, 1u);
  ASSERT_EQ(n_overlap, 1u);
  strip_pair_collection_types::host standard(n_standard, &mr);
  strip_pair_collection_types::host overlap(n_overlap, &mr);
  for (unsigned int i = 0; i < 4; ++i) {
    device::find_strip_pairs(
        i, meas_view, surfaces_view, rules_view, point3{},
        vecmem::get_data(standard_offsets), vecmem::get_data(overlap_offsets),
        vecmem::get_data(standard), vecmem::get_data(overlap));
  }
  EXPECT_EQ(standard[0].measurement_index_2, 1u);
  EXPECT_EQ(overlap[0].measurement_index_2, 2u);
  // Deliberately reserve no space for measurement zero. Find must not write
  // into the following measurement's reserved slot, even if it finds a pair.
  const strip_pair sentinel{std::numeric_limits<unsigned int>::max(),
                            std::numeric_limits<unsigned int>::max(), 0.f, 0.f};
  strip_pair_collection_types::host guarded_standard(1u, sentinel, &mr);
  strip_pair_collection_types::host guarded_overlap(1u, sentinel, &mr);
  vecmem::vector<unsigned int> underestimated_offsets{&mr};
  underestimated_offsets.assign({0u, 1u, 1u});
  device::find_strip_pairs(0u, meas_view, surfaces_view, rules_view, point3{},
                           vecmem::get_data(underestimated_offsets),
                           vecmem::get_data(underestimated_offsets),
                           vecmem::get_data(guarded_standard),
                           vecmem::get_data(guarded_overlap));
  EXPECT_EQ(guarded_standard[0].measurement_index_1,
            sentinel.measurement_index_1);
  EXPECT_EQ(guarded_overlap[0].measurement_index_1,
            sentinel.measurement_index_1);
  for (const auto* pairs : {&standard, &overlap}) {
    strip_pair_collection_types::host pairs_with_sentinel{&mr};
    pairs_with_sentinel.push_back(pairs->at(0));
    pairs_with_sentinel.push_back({std::numeric_limits<unsigned int>::max(),
                                   std::numeric_limits<unsigned int>::max(),
                                   0.f, 0.f});
    edm::spacepoint_collection::buffer candidates(2u, mr);
    vecmem::vector<unsigned int> accepted(2u, 0u, &mr);
    copy.setup(candidates)->wait();
    // Process a valid pair, an invalid sentinel, and an out-of-range lane.
    for (unsigned int i = 0; i < 3; ++i) {
      device::form_strip_spacepoints_from_pairs(
          i, meas_view, vecmem::get_data(pairs_with_sentinel), surfaces_view,
          point3{}, vecmem::get_data(accepted), candidates);
    }
    EXPECT_EQ(accepted[0], 1u);
    EXPECT_EQ(accepted[1], 0u);
    std::partial_sum(accepted.begin(), accepted.end(), accepted.begin());
    edm::spacepoint_collection::buffer output(accepted.back(), mr);
    copy.setup(output)->wait();
    const edm::spacepoint_collection::const_view candidates_view(candidates);
    // Gather in reverse execution order to verify scanned destination indices.
    for (unsigned int i = 3u; i > 0u; --i) {
      device::gather_strip_spacepoints(i - 1u, candidates_view,
                                       vecmem::get_data(accepted), output);
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
