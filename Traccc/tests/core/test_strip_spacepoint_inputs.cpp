// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <gtest/gtest.h>
#include <vecmem/memory/host_memory_resource.hpp>

#include "traccc/seeding/detail/spacepoint_formation.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"
#include "traccc/seeding/detail/strip_pairing.hpp"

namespace {
using namespace traccc;
using measurements_t = edm::measurement_collection::host;

detray::geometry::identifier barcode(unsigned int index,
                                     unsigned int volume = 0u) {
  detray::geometry::identifier result{0u};
  result.set_index(index);
  result.set_volume(volume);
  return result;
}

void add_measurement(measurements_t& measurements, unsigned int index,
                     float local0 = 0.f, float local1 = 0.f,
                     unsigned int volume = 0u, unsigned int dimensions = 1u) {
  measurements.push_back({{local0, local1},
                          {0.f, 0.f},
                          dimensions,
                          0.f,
                          0.f,
                          0u,
                          barcode(index, volume),
                          {1u, 0u},
                          0u});
}

strip_measurement_surface_info linear_info(std::uint64_t link) {
  strip_measurement_surface_info info{};
  info.surface_link = link;
  info.origin = {100.f, 0.f, 0.f};
  info.local_u = {0.f, 0.f, 1.f};
  info.local_v = {0.f, 1.f, 0.f};
  info.linear = {1u, 10u, 1.f, 20.f, 0.f};
  info.index_mapping = {0u, 10u, 1.f, 1.f, 5.f, false};
  return info;
}

struct pair_collector {
  std::vector<unsigned int> indices;
  void operator()(unsigned int index, const strip_pairing_rule&) {
    indices.push_back(index);
  }
};

TEST(strip_spacepoint_inputs, index_rounding_and_invalid_input) {
  vecmem::host_memory_resource mr;
  measurements_t measurements{mr};
  add_measurement(measurements, 0u, -3.f);
  const auto measurement = measurements.at(0u);
  // Odd linear strip counts truncate the half-integer before range checks.
  strip_index_mapping mapping{0u, 5u, 1.f, 1.f, 2.5f, false};
  EXPECT_EQ(details::measured_strip_index(measurement, mapping), 0);
  mapping.check_before_truncation = true;
  EXPECT_EQ(details::measured_strip_index(measurement, mapping), -1);

  measurements.at(0u).local_position()[1] = 0.f;
  mapping = {1u, 10u, 0.01f, -1.f, 4.5f, true};
  EXPECT_EQ(details::measured_strip_index(measurement, mapping), 4);
  mapping.pitch = 0.f;
  EXPECT_EQ(details::measured_strip_index(measurement, mapping), -1);
  mapping.pitch = 0.01f;
  measurements.at(0u).local_position()[1] =
      std::numeric_limits<scalar>::quiet_NaN();
  EXPECT_EQ(details::measured_strip_index(measurements.at(0u), mapping), -1);
}

TEST(strip_spacepoint_inputs, linear_frame_and_event_beamspot) {
  vecmem::host_memory_resource mr;
  measurements_t measurements{mr};
  add_measurement(measurements, 0u);
  auto info = linear_info(barcode(0u).value());
  const auto material =
      details::make_strip_material(measurements.at(0u), info, point3{});
  ASSERT_EQ(material.valid, 1u);
  EXPECT_EQ(material.center[0], scalar{100});
  EXPECT_EQ(material.center[1], scalar{0.5f});
  EXPECT_EQ(material.direction[2], scalar{-20});
  EXPECT_EQ(material.half_length, scalar{10});
  const auto shifted = details::make_strip_material(measurements.at(0u), info,
                                                    point3{1.f, 2.f, 3.f});
  EXPECT_EQ(shifted.center[0], material.center[0]);
  EXPECT_EQ(shifted.trajectory[0], material.trajectory[0] - scalar{2});
  EXPECT_EQ(shifted.trajectory[1], material.trajectory[1] - scalar{4});
  EXPECT_EQ(shifted.trajectory[2], material.trajectory[2] - scalar{6});
  // A different frame orientation requires no barrel/endcap classification.
  info.local_u = {1.f, 0.f, 0.f};
  info.local_v = {0.f, 0.f, 1.f};
  const auto rotated =
      details::make_strip_material(measurements.at(0u), info, point3{});
  EXPECT_EQ(rotated.direction[0], scalar{-20});
  EXPECT_EQ(rotated.direction[2], scalar{0});
}

TEST(strip_spacepoint_inputs, radial_geometry_without_region_flag) {
  radial_strip_geometry geometry{};
  geometry.n_strips = 1u;
  geometry.angular_pitch = 0.01f;
  geometry.min_radius = 10.f;
  geometry.max_radius = 20.f;
  const auto first =
      details::radial_strip_position_at_radius(0, 10.f, geometry);
  const auto second =
      details::radial_strip_position_at_radius(0, 20.f, geometry);
  EXPECT_FLOAT_EQ(first[0], 10.f);
  EXPECT_FLOAT_EQ(first[1], 0.f);
  EXPECT_FLOAT_EQ(second[0], 20.f);
  geometry.local_frame = strip_local_frame::polar;
  const auto polar =
      details::radial_strip_position_at_radius(0, 10.f, geometry);
  EXPECT_FLOAT_EQ(polar[0], 10.f);
  EXPECT_FLOAT_EQ(polar[1], 0.f);
}

TEST(strip_spacepoint_inputs, signed_difference_and_inclusive_windows) {
  vecmem::host_memory_resource mr;
  measurements_t measurements{mr};
  add_measurement(measurements, 0u, 1.f);
  add_measurement(measurements, 1u, 3.f);
  auto reference = linear_info(barcode(0u).value());
  auto candidate = linear_info(barcode(1u).value());
  strip_pairing_rule rule{};
  rule.reference_surface_link = reference.surface_link;
  rule.candidate_surface_link = candidate.surface_link;
  rule.difference_scale = -1.f;
  rule.difference_min = rule.difference_max = -2.f;
  EXPECT_TRUE(details::match_strip_pair(
      measurements.at(0u), measurements.at(1u), reference, candidate, rule));
  rule.difference_min = rule.difference_max = 2.f;
  EXPECT_FALSE(details::match_strip_pair(
      measurements.at(0u), measurements.at(1u), reference, candidate, rule));
  rule.mode = strip_pairing_mode::windows;
  rule.reference_min = rule.reference_max = 1.f;
  rule.candidate_min = rule.candidate_max = 3.f;
  EXPECT_TRUE(details::match_strip_pair(
      measurements.at(0u), measurements.at(1u), reference, candidate, rule));
  rule.reference_coordinate = rule.candidate_coordinate =
      strip_pairing_coordinate::strip_index;
  rule.reference_min = rule.reference_max = 6.f;
  rule.candidate_min = rule.candidate_max = 8.f;
  EXPECT_TRUE(details::match_strip_pair(
      measurements.at(0u), measurements.at(1u), reference, candidate, rule));
}

TEST(strip_spacepoint_inputs, arbitrary_rule_count_and_barcode_order) {
  vecmem::host_memory_resource mr;
  measurements_t measurements{mr};
  std::vector<strip_measurement_surface_info> surfaces;
  std::vector<strip_pairing_rule> rules;
  // The reference has the smallest index but largest full barcode.
  add_measurement(measurements, 0u, 0.f, 0.f, 9u);
  const auto reference_link = barcode(0u, 9u).value();
  surfaces.push_back(linear_info(reference_link));
  for (unsigned int i = 1u; i <= 7u; ++i) {
    add_measurement(measurements, i);
    surfaces.push_back(linear_info(barcode(i).value()));
    strip_pairing_rule rule{};
    rule.reference_surface_link = reference_link;
    rule.candidate_surface_link = barcode(i).value();
    rules.push_back(rule);
  }
  // A measurement with the same surface index but a different full barcode
  // must not match the rule. Two-dimensional measurements must not match.
  add_measurement(measurements, 7u, 0.f, 0.f, 2u);
  add_measurement(measurements, 7u, 0.f, 0.f, 0u, 2u);
  std::sort(surfaces.begin(), surfaces.end(), [](const auto& a, const auto& b) {
    return a.surface_link < b.surface_link;
  });
  std::sort(rules.begin(), rules.end(), [](const auto& a, const auto& b) {
    return a.candidate_surface_link < b.candidate_surface_link;
  });
  pair_collector count_pass, write_pass;
  const point3 beam_spot{1.f, 2.f, 3.f};
  details::visit_strip_pairs(0u, measurements, surfaces, rules, beam_spot,
                             count_pass);
  details::visit_strip_pairs(0u, measurements, surfaces, rules, beam_spot,
                             write_pass);
  EXPECT_EQ(count_pass.indices.size(), 7u);
  EXPECT_EQ(count_pass.indices, write_pass.indices);
  pair_collector reverse;
  details::visit_strip_pairs(1u, measurements, surfaces, rules, beam_spot,
                             reverse);
  EXPECT_TRUE(reverse.indices.empty());
  // Both passes reject an invalid candidate geometry identically.
  surfaces.front().index_mapping.pitch = 0.f;
  pair_collector invalid;
  details::visit_strip_pairs(0u, measurements, surfaces, rules, beam_spot,
                             invalid);
  EXPECT_EQ(invalid.indices.size(), 6u);
}

TEST(strip_spacepoint_inputs, intersection_and_degenerate_pair) {
  point3 position{};
  // Two strips cross the beam ray at (100,0,0).
  EXPECT_TRUE(details::make_strip_spacepoint(
      position, point3{100.f, 0.f, 0.f}, vector3{0.f, 0.f, 20.f},
      vector3{0.f, 20.f, 0.f}, vector3{200.f, 0.f, 0.f},
      vector3{200.f, 0.f, 0.f}, vector3{0.f, 4000.f, 0.f},
      vector3{0.f, 0.f, -4000.f}, 10.f, 10.f, 0.f, 0.01f));
  EXPECT_EQ(position[0], scalar{100});
  EXPECT_EQ(position[1], scalar{0});
  EXPECT_EQ(position[2], scalar{0});
  // Orthogonal strips have cs == 0, but require no endpoint correction.
  // A nonzero gap allowance must not reject this valid intersection.
  EXPECT_TRUE(details::make_strip_spacepoint(
      position, point3{100.f, 0.f, 0.f}, vector3{0.f, 0.f, 20.f},
      vector3{0.f, 20.f, 0.f}, vector3{200.f, 0.f, 0.f},
      vector3{200.f, 0.f, 0.f}, vector3{0.f, 4000.f, 0.f},
      vector3{0.f, 0.f, -4000.f}, 10.f, 10.f, 1.f, 0.01f));
  EXPECT_FALSE(details::make_strip_spacepoint(
      position, point3{100.f, 0.f, 0.f}, vector3{0.f, 0.f, 20.f},
      vector3{0.f, 0.f, 20.f}, vector3{200.f, 0.f, 0.f},
      vector3{200.f, 0.f, 0.f}, vector3{0.f, 4000.f, 0.f},
      vector3{0.f, 4000.f, 0.f}, 10.f, 10.f, 0.f, 0.01f));
}
}  // namespace
