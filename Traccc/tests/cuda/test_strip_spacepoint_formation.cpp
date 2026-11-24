// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <array>
#include <vector>

#include <gtest/gtest.h>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

#include "traccc/cuda/seeding/silicon_strip_spacepoint_formation_algorithm.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"

namespace {
using namespace traccc;

class CUDAStripSpacepointFormation : public ::testing::Test {
 protected:
  vecmem::cuda::device_memory_resource device_mr;
  vecmem::cuda::host_memory_resource host_mr;
  vecmem::cuda::stream_wrapper owning_stream;
  cuda::stream_wrapper stream{owning_stream.stream()};
  vecmem::cuda::async_copy copy{stream.cudaStream()};
  memory_resource mr{device_mr, &host_mr};
  edm::measurement_collection::host measurements{host_mr};
  strip_measurement_surface_info_collection_types::host surfaces{&host_mr};
  strip_pairing_rule_collection_types::host rules{&host_mr};

  struct point {
    std::array<float, 3> position;
    float variance_r, variance_z;
    unsigned int first, second;
  };
  struct result {
    std::vector<point> standard, overlap;
    std::vector<strip_spacepoint_shape> shapes, overlap_shapes;
  };

  void SetUp() override {
    // Adapter geometry: orthogonal strips cross at (100, 0, 0).
    // Use synthetic identifiers; Strip formation consumes the adapter geometry.
    for (unsigned int i = 0; i < 3; ++i) {
      detray::geometry::identifier link{0u};
      link.set_index(i);
      measurements.push_back({{0.f, 0.f},
                              {0.f, 0.f},
                              1u,
                              0.f,
                              0.f,
                              0u,
                              link,
                              {1u, 0u},
                              i,
                              {0.f, 0.f}});
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
        rule.category = i == 1 ? strip_pair_category::standard
                               : strip_pair_category::overlap;
        rules.push_back(rule);
      }
    }
  }

  result run(const point3& beam_spot = {}, bool retain_shapes = true) {
    auto device_measurements =
        copy.to(vecmem::get_data(measurements), device_mr, &host_mr,
                vecmem::copy::type::host_to_device);
    auto device_surfaces = copy.to(vecmem::get_data(surfaces), device_mr,
                                   vecmem::copy::type::host_to_device);
    auto device_rules = copy.to(vecmem::get_data(rules), device_mr,
                                vecmem::copy::type::host_to_device);
    cuda::silicon_strip_spacepoint_formation_algorithm algorithm{mr, copy,
                                                                 stream};
    auto output = retain_shapes
                      ? algorithm(device_measurements, device_surfaces,
                                  device_rules, beam_spot, true)
                      : algorithm(device_measurements, device_surfaces,
                                  device_rules, beam_spot);
    auto read = [&](const auto& buffer) {
      std::vector<point> points;
      const auto size = copy.get_size(buffer);
      if (size != 0u) {
        auto host_buffer = copy.to(buffer, host_mr, nullptr,
                                   vecmem::copy::type::device_to_host);
        stream.synchronize();
        edm::spacepoint_collection::const_device host_points{host_buffer};
        for (unsigned int i = 0; i < size; ++i) {
          const auto sp = host_points.at(i);
          points.push_back({sp.global(), sp.radius_variance(), sp.z_variance(),
                            sp.measurement_index_1(),
                            sp.measurement_index_2()});
        }
      }
      return points;
    };
    result answer{};
    answer.standard = read(output.spacepoints);
    answer.overlap = read(output.overlap_spacepoints);
    const auto read_shapes = [&](const auto& buffer) {
      std::vector<strip_spacepoint_shape> shapes;
      const auto size = copy.get_size(buffer);
      if (size != 0u) {
        auto host =
            copy.to(buffer, host_mr, vecmem::copy::type::device_to_host);
        stream.synchronize();
        const strip_spacepoint_shape_collection_types::const_device values(
            host);
        for (unsigned int i = 0; i < size; ++i)
          shapes.push_back(values.at(i));
      }
      return shapes;
    };
    answer.shapes = read_shapes(output.shapes);
    answer.overlap_shapes = read_shapes(output.overlap_shapes);
    EXPECT_EQ(answer.shapes.size(),
              retain_shapes ? answer.standard.size() : 0u);
    EXPECT_EQ(answer.overlap_shapes.size(),
              retain_shapes ? answer.overlap.size() : 0u);
    const auto check_shapes = [&](const auto& points, const auto& shapes) {
      for (std::size_t i = 0; i < shapes.size(); ++i) {
        const auto first = measurements.at(points[i].first);
        const auto second = measurements.at(points[i].second);
        strip_measurement_surface_info a{}, b{};
        ASSERT_TRUE(details::find_strip_surface_info(
            surfaces, first.surface_link().value(), a));
        ASSERT_TRUE(details::find_strip_surface_info(
            surfaces, second.surface_link().value(), b));
        const auto top = details::make_strip_material(first, a, beam_spot);
        const auto bottom = details::make_strip_material(second, b, beam_spot);
        EXPECT_NEAR(shapes[i].top_half_length, top.half_length, 1e-5);
        EXPECT_NEAR(shapes[i].bottom_half_length, bottom.half_length, 1e-5);
        for (unsigned int k = 0; k < 3u; ++k) {
          EXPECT_NEAR(shapes[i].top_direction[k],
                      -top.direction[k] / (2 * top.half_length), 1e-5);
          EXPECT_NEAR(shapes[i].bottom_direction[k],
                      -bottom.direction[k] / (2 * bottom.half_length), 1e-5);
          EXPECT_NEAR(shapes[i].top_center[k], top.center[k], 1e-5);
          EXPECT_NEAR(shapes[i].center_distance[k],
                      top.center[k] - bottom.center[k], 1e-5);
        }
      }
    };
    check_shapes(answer.standard, answer.shapes);
    check_shapes(answer.overlap, answer.overlap_shapes);
    stream.synchronize();
    return answer;
  }

  void check_point(const point& p, unsigned int first, unsigned int second) {
    EXPECT_FLOAT_EQ(p.position[0], 100.f);
    EXPECT_FLOAT_EQ(p.position[1], 0.f);
    EXPECT_FLOAT_EQ(p.position[2], 0.f);
    EXPECT_FLOAT_EQ(p.variance_r, 0.1f);
    EXPECT_FLOAT_EQ(p.variance_z, 0.2f);
    EXPECT_EQ(p.first, first);
    EXPECT_EQ(p.second, second);
  }
};

TEST_F(CUDAStripSpacepointFormation, standard_and_overlap) {
  const auto output = run();
  ASSERT_EQ(output.standard.size(), 1u);
  ASSERT_EQ(output.overlap.size(), 1u);
  check_point(output.standard[0], 0u, 1u);
  check_point(output.overlap[0], 0u, 2u);
}

TEST_F(CUDAStripSpacepointFormation, optional_shapes_preserve_spacepoints) {
  const auto legacy = run({}, false);
  const auto retained = run({}, true);
  ASSERT_EQ(legacy.standard.size(), retained.standard.size());
  ASSERT_EQ(legacy.overlap.size(), retained.overlap.size());
  EXPECT_TRUE(legacy.shapes.empty());
  EXPECT_TRUE(legacy.overlap_shapes.empty());
  for (const auto& values : {std::pair{&legacy.standard, &retained.standard},
                             std::pair{&legacy.overlap, &retained.overlap}}) {
    for (std::size_t i = 0; i < values.first->size(); ++i) {
      EXPECT_EQ(values.first->at(i).position, values.second->at(i).position);
      EXPECT_EQ(values.first->at(i).first, values.second->at(i).first);
      EXPECT_EQ(values.first->at(i).second, values.second->at(i).second);
    }
  }
}

TEST_F(CUDAStripSpacepointFormation, radial_shapes_follow_endpoint_convention) {
  for (auto& surface : surfaces) {
    surface.geometry_model = strip_geometry_model::radial;
    surface.radial.n_strips = 1u;
    surface.radial.angular_pitch = 0.01f;
    surface.radial.min_radius = 10.f;
    surface.radial.max_radius = 20.f;
    for (unsigned int k = 0; k < 3u; ++k)
      surface.origin[k] = (k == 0u ? 100.f : 0.f) - 15.f * surface.local_u[k];
  }
  const auto output = run();
  ASSERT_EQ(output.standard.size(), 1u);
  ASSERT_EQ(output.overlap.size(), 1u);
  check_point(output.standard[0], 0u, 1u);
  EXPECT_NEAR(output.shapes[0].top_half_length, 5.f, 1e-5);
  EXPECT_NEAR(output.shapes[0].bottom_half_length, 5.f, 1e-5);
  EXPECT_NEAR(output.shapes[0].top_direction[2], 1.f, 1e-5);
  EXPECT_NEAR(output.shapes[0].bottom_direction[1], 1.f, 1e-5);
}

TEST_F(CUDAStripSpacepointFormation, standard_only) {
  rules.resize(1);
  const auto output = run();
  ASSERT_EQ(output.standard.size(), 1u);
  EXPECT_TRUE(output.overlap.empty());
  check_point(output.standard[0], 0u, 1u);
}

TEST_F(CUDAStripSpacepointFormation, empty_measurements) {
  measurements.resize(0);
  const auto output = run();
  EXPECT_TRUE(output.standard.empty());
  EXPECT_TRUE(output.overlap.empty());
}

TEST_F(CUDAStripSpacepointFormation, no_pairing_rules) {
  rules.clear();
  const auto output = run();
  EXPECT_TRUE(output.standard.empty());
  EXPECT_TRUE(output.overlap.empty());
}

TEST_F(CUDAStripSpacepointFormation, rejected_pair_window) {
  for (auto& rule : rules) {
    rule.difference_min = rule.difference_max = 1.f;
  }
  const auto output = run();
  EXPECT_TRUE(output.standard.empty());
  EXPECT_TRUE(output.overlap.empty());
}

TEST_F(CUDAStripSpacepointFormation, parallel_strips) {
  // Pairs pass the selection, but have no valid spacepoint intersection.
  for (auto& surface : surfaces) {
    surface.local_u = {0.f, 0.f, 1.f};
    surface.local_v = {0.f, 1.f, 0.f};
  }
  const auto output = run();
  EXPECT_TRUE(output.standard.empty());
  EXPECT_TRUE(output.overlap.empty());
}

TEST_F(CUDAStripSpacepointFormation, pixel_measurement_is_ignored) {
  measurements.at(1).dimensions() = 2u;
  const auto output = run();
  EXPECT_TRUE(output.standard.empty());
  ASSERT_EQ(output.overlap.size(), 1u);
  check_point(output.overlap[0], 0u, 2u);
}

TEST_F(CUDAStripSpacepointFormation, shifted_beam_spot) {
  const auto output = run(point3{1.f, 2.f, 3.f});
  ASSERT_EQ(output.standard.size(), 1u);
  ASSERT_EQ(output.overlap.size(), 1u);
  check_point(output.standard[0], 0u, 1u);
  check_point(output.overlap[0], 0u, 2u);
}

TEST_F(CUDAStripSpacepointFormation, multiple_blocks_and_repeated_calls) {
  constexpr unsigned int n_reference = 257u;
  measurements.resize(0);
  for (unsigned int i = 0; i < n_reference + 2u; ++i) {
    const unsigned int surface = i < n_reference ? 0u : i - n_reference + 1u;
    measurements.push_back(
        {{0.f, 0.f},
         {0.f, 0.f},
         1u,
         0.f,
         0.f,
         0u,
         detray::geometry::identifier{surfaces[surface].surface_link},
         {1u, 0u},
         i,
         {0.f, 0.f}});
  }
  for (unsigned int repeat = 0; repeat < 3; ++repeat) {
    auto output = run();
    ASSERT_EQ(output.standard.size(), n_reference);
    ASSERT_EQ(output.overlap.size(), n_reference);
    // Check the native output order; do not sort away scheduling differences.
    for (unsigned int i = 0; i < n_reference; ++i) {
      check_point(output.standard[i], i, n_reference);
      check_point(output.overlap[i], i, n_reference + 1u);
    }
  }
}
}  // namespace
