// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/builders/detail/bounding_box_creator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/common/build_nested_detector.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace detray;

namespace {

using algebra_t = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

/// @returns a printable string for the point @p p
std::string to_string(const point3 &p) {
  return "[" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " +
         std::to_string(p[2]) + "]";
}

/// A box of the nested test detector
struct box {
  /// Volume index of the box
  dindex index;
  /// Center of the box
  point3 center;
  /// Half lengths of the box
  std::array<scalar, 3> half_lengths;

  /// @returns the signed distance of @p p to the box boundary
  /// (positive: inside)
  scalar inner_distance(const point3 &p) const {
    scalar dist{std::numeric_limits<scalar>::max()};
    for (unsigned int i = 0u; i < 3u; ++i) {
      dist = std::min(dist, half_lengths[i] - std::abs(p[i] - center[i]));
    }
    return dist;
  }
};

/// The boxes of the nested test detector (no inset), mothers before
/// daughters
const std::vector<box> nested_boxes{
    {0u, {0.f, 0.f, 0.f}, {200.f, 200.f, 200.f}},
    {1u, {0.f, 0.f, -100.f}, {60.f, 60.f, 40.f}},
    {2u, {0.f, 0.f, -130.f}, {20.f, 20.f, 10.f}},
    {3u, {0.f, 0.f, 160.f}, {60.f, 60.f, 40.f}},
    {4u, {0.f, 130.f, 0.f}, {40.f, 40.f, 30.f}},
    {5u, {0.f, 130.f, 60.f}, {40.f, 40.f, 30.f}}};

/// Count the volumes that contain a point
template <typename detector_t>
void count_containing_volumes(const detector_t &det, const point3 &p,
                              unsigned int &n_inside, unsigned int &n_boundary,
                              unsigned int &n_ambiguous) {
  n_inside = 0u;
  n_boundary = 0u;
  n_ambiguous = 0u;
  for (const auto &vol_desc : det.volumes()) {
    switch (tracking_volume{det, vol_desc}.is_inside(p)) {
      case volume_containment::e_inside:
        ++n_inside;
        break;
      case volume_containment::e_on_boundary:
        ++n_boundary;
        break;
      case volume_containment::e_ambiguous:
        ++n_ambiguous;
        break;
      default:
        break;
    }
  }
}

/// Every point inside the detector must be inside exactly one volume (flat
/// detectors fill the space without overlaps)
template <typename detector_t>
void test_flat_detector_containment(const detector_t &det,
                                    const unsigned int n_points) {
  // Extent of the detector
  point3 det_min{std::numeric_limits<scalar>::max(),
                 std::numeric_limits<scalar>::max(),
                 std::numeric_limits<scalar>::max()};
  point3 det_max{-det_min[0], -det_min[1], -det_min[2]};
  for (const auto &vol_desc : det.volumes()) {
    const auto aabb{detail::volume_bounding_box(det, vol_desc.index())};
    const point3 v_min{aabb.template loc_min<point3>()};
    const point3 v_max{aabb.template loc_max<point3>()};
    for (unsigned int i = 0u; i < 3u; ++i) {
      det_min[i] = std::min(det_min[i], v_min[i]);
      det_max[i] = std::max(det_max[i], v_max[i]);
    }
  }

  // The detectors are cylinders around the z-axis
  const scalar r_max{0.99f * det_max[0]};
  const scalar z_min{0.99f * det_min[2]};
  const scalar z_max{0.99f * det_max[2]};

  std::mt19937 gen{42u};
  std::uniform_real_distribution<scalar> dist_r{0.f, 1.f};
  std::uniform_real_distribution<scalar> dist_phi{-constant<scalar>::pi,
                                                  constant<scalar>::pi};
  std::uniform_real_distribution<scalar> dist_z{z_min, z_max};

  unsigned int n_ambiguous_total{0u};
  for (unsigned int i = 0u; i < n_points; ++i) {
    // Uniform in the cylinder volume
    const scalar r{r_max * std::sqrt(dist_r(gen))};
    const scalar phi{dist_phi(gen)};
    const point3 p{r * std::cos(phi), r * std::sin(phi), dist_z(gen)};

    unsigned int n_inside{0u};
    unsigned int n_boundary{0u};
    unsigned int n_ambiguous{0u};
    count_containing_volumes(det, p, n_inside, n_boundary, n_ambiguous);

    n_ambiguous_total += n_ambiguous;
    if (n_ambiguous == 0u && n_boundary == 0u) {
      EXPECT_EQ(n_inside, 1u) << "at point " << to_string(p);
    }
    // Never inside of more than one volume
    EXPECT_LE(n_inside, 1u) << "at point " << to_string(p);
  }
  // Cylinder volumes are never ambiguous (no ray test)
  EXPECT_EQ(n_ambiguous_total, 0u);
}

/// Build a detector with a single, rotated cuboid volume
///
/// @param mr memory resource for the detector containers
/// @param center the center of the box
/// @param axes the (orthonormal) axes of the box
/// @param half_lengths the half lengths of the box along its axes
/// @param scale scale factor for the size of the first portal
/// @param all_parallel make all portals parallel (invalid cuboid)
auto build_rotated_box(vecmem::memory_resource &mr, const point3 &center,
                       const std::array<test::vector3, 3> &axes,
                       const std::array<scalar, 3> &half_lengths,
                       const scalar scale = 1.f,
                       const bool all_parallel = false) {
  using builder_t = detector_builder<nested_volumes_metadata<algebra_t>>;
  using detector_t = typename builder_t::detector_type;
  using transform3_t = typename detector_t::transform3_type;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  constexpr nav_link_t end_of_world{detail::invalid_value<nav_link_t>()};

  builder_t det_builder{};
  auto *v_builder = det_builder.new_volume(volume_id::e_cuboid);
  v_builder->add_volume_placement(center, axes[0], axes[2]);

  // Add the faces in an unsorted order: (axis, side)
  auto pt_factory =
      std::make_shared<surface_factory<detector_t, rectangle2D>>();
  const std::array<std::pair<unsigned int, scalar>, 6> faces{
      {{2u, 1.f}, {0u, -1.f}, {2u, -1.f}, {1u, 1.f}, {0u, 1.f}, {1u, -1.f}}};
  for (const auto &[i, face] : detray::views::enumerate(faces)) {
    const unsigned int k{all_parallel ? 2u : face.first};
    const scalar side{face.second};
    const unsigned int kx{(k + 1u) % 3u};
    const unsigned int ky{(k + 2u) % 3u};

    // Local x-axis along the next box axis, normal along the box axis
    const test::vector3 n{side * axes[k]};
    const scalar pos{all_parallel ? static_cast<scalar>(i) : half_lengths[k]};
    const transform3_t trf{center + pos * side * axes[k], n, axes[kx]};

    const scalar f{i == 0u ? scale : 1.f};
    const scalar hx{f * half_lengths[kx]};
    const scalar hy{f * half_lengths[ky]};
    pt_factory->push_back({surface_id::e_portal, trf, end_of_world,
                           std::vector<scalar>{-hx, -hy, hx, hy}});
  }
  v_builder->add_surfaces(pt_factory);

  return det_builder.build(mr, volume_builder_options{});
}

}  // anonymous namespace

/// Unittest: Volume containment for the boxes of the nested detector
GTEST_TEST(detray_geometry, volume_containment_nested_boxes) {
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = build_nested_detector<algebra_t>(host_mr);

  const auto world = tracking_volume{det, 0u};
  const auto box_a = tracking_volume{det, 1u};

  // Inside (daughter volumes are part of the mother volume)
  EXPECT_EQ(world.is_inside({1.f, 2.f, 3.f}), volume_containment::e_inside);
  EXPECT_EQ(world.is_inside({1.f, 2.f, -100.f}), volume_containment::e_inside);
  EXPECT_EQ(box_a.is_inside({1.f, 2.f, -100.f}), volume_containment::e_inside);

  // Outside
  EXPECT_EQ(world.is_inside({1.f, 2.f, 250.f}), volume_containment::e_outside);
  EXPECT_EQ(box_a.is_inside({1.f, 2.f, 3.f}), volume_containment::e_outside);
  EXPECT_EQ(box_a.is_inside({70.f, 0.f, -100.f}),
            volume_containment::e_outside);

  // On a face, an edge and a corner
  EXPECT_EQ(world.is_inside({200.f, 3.f, 5.f}),
            volume_containment::e_on_boundary);
  EXPECT_EQ(box_a.is_inside({3.f, 5.f, -60.f}),
            volume_containment::e_on_boundary);
  EXPECT_EQ(box_a.is_inside({60.f, 60.f, -100.f}),
            volume_containment::e_on_boundary);
  EXPECT_EQ(box_a.is_inside({60.f, 60.f, -60.f}),
            volume_containment::e_on_boundary);

  // Close to a face, but outside of the tolerance
  EXPECT_EQ(box_a.is_inside({3.f, 5.f, -60.01f}), volume_containment::e_inside);
  EXPECT_EQ(box_a.is_inside({3.f, 5.f, -59.99f}),
            volume_containment::e_outside);
}

/// Unittest: Find the volume of a point in the nested detector
GTEST_TEST(detray_geometry, volume_lookup_nested_detector) {
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = build_nested_detector<algebra_t>(host_mr);

  // Do not test points close to a box boundary
  constexpr scalar margin{0.01f};

  unsigned int n_tested{0u};
  std::array<unsigned int, 6> n_per_volume{};
  for (scalar x = -213.7f; x < 214.f; x += 7.3f) {
    for (scalar y = -213.1f; y < 214.f; y += 7.9f) {
      for (scalar z = -212.9f; z < 214.f; z += 6.7f) {
        const point3 p{x, y, z};

        // Truth: The innermost box that contains the point
        dindex exp_volume{dindex_invalid};
        bool is_close{false};
        for (const auto &b : nested_boxes) {
          const scalar dist{b.inner_distance(p)};
          is_close |= std::abs(dist) < margin;
          if (dist > 0.f) {
            exp_volume = b.index;
          }
        }
        if (is_close) {
          continue;
        }

        ++n_tested;
        const dindex vol_idx{det.volume_index(p)};
        ASSERT_EQ(vol_idx, exp_volume) << "at point " << to_string(p);

        if (vol_idx != dindex_invalid) {
          ++n_per_volume.at(vol_idx);
          EXPECT_EQ(det.volume(p).index(), vol_idx);
        }
      }
    }
  }

  // Every volume was tested
  EXPECT_GT(n_tested, 100000u);
  for (const unsigned int n : n_per_volume) {
    EXPECT_GT(n, 10u);
  }
}

/// Unittest: Volume containment in the toy detector (cylinder volumes)
GTEST_TEST(detray_geometry, volume_containment_toy_detector) {
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = build_toy_detector<algebra_t>(host_mr);

  test_flat_detector_containment(det, 2000u);

  // The top-level lookup without a volume accelerator finds the same volume
  const point3 p{12.3f, -4.5f, 67.8f};
  const dindex vol_idx{det.volume_index(p)};
  ASSERT_NE(vol_idx, dindex_invalid);
  EXPECT_EQ(tracking_volume(det, vol_idx).is_inside(p),
            volume_containment::e_inside);
}

/// Unittest: Volume containment in the wire chamber (cylinder volumes)
GTEST_TEST(detray_geometry, volume_containment_wire_chamber) {
  vecmem::host_memory_resource host_mr;
  wire_chamber_config<scalar> wire_cfg{};
  const auto [det, names] = build_wire_chamber<algebra_t>(host_mr, wire_cfg);

  test_flat_detector_containment(det, 2000u);
}

/// Unittest: Inside check for a rotated cuboid volume
GTEST_TEST(detray_geometry, volume_containment_rotated_box) {
  using vector3 = test::vector3;

  vecmem::host_memory_resource host_mr;

  // Box axes: Rotated around the z-axis and the (new) x-axis
  const scalar a{0.5f};
  const scalar b{0.3f};
  const vector3 e0{std::cos(a), std::sin(a), 0.f};
  const vector3 e1_r{-std::sin(a), std::cos(a), 0.f};
  const vector3 e1{std::cos(b) * e1_r + std::sin(b) * vector3{0.f, 0.f, 1.f}};
  const vector3 e2{vector::cross(e0, e1)};
  const std::array<vector3, 3> axes{e0, e1, e2};
  const point3 center{10.f, -20.f, 30.f};
  const std::array<scalar, 3> h{10.f, 20.f, 30.f};

  const auto det = build_rotated_box(host_mr, center, axes, h);
  EXPECT_NO_THROW(detail::check_consistency(det));

  const auto vol = tracking_volume{det, 0u};
  ASSERT_EQ(vol.n_portals(), 6u);

  // Points on a face and on a corner
  auto to_global = [&](const scalar l0, const scalar l1, const scalar l2) {
    return center + l0 * e0 + l1 * e1 + l2 * e2;
  };
  EXPECT_EQ(vol.is_inside(to_global(h[0], 0.f, 0.f)),
            volume_containment::e_on_boundary);
  EXPECT_EQ(vol.is_inside(to_global(-h[0], -h[1], h[2])),
            volume_containment::e_on_boundary);
  EXPECT_EQ(vol.is_inside(to_global(0.f, 0.f, 0.f)),
            volume_containment::e_inside);

  // Random points around the box
  std::mt19937 gen(42u);
  std::uniform_real_distribution<scalar> dist(-1.5f, 1.5f);
  constexpr scalar margin{0.01f};
  for (unsigned int i = 0u; i < 10000u; ++i) {
    const std::array<scalar, 3> l{dist(gen) * h[0], dist(gen) * h[1],
                                  dist(gen) * h[2]};
    scalar inner_dist{std::numeric_limits<scalar>::max()};
    for (unsigned int k = 0u; k < 3u; ++k) {
      inner_dist = std::min(inner_dist, h[k] - std::abs(l[k]));
    }
    if (std::abs(inner_dist) < margin) {
      continue;
    }
    const auto exp_vc{inner_dist > 0.f ? volume_containment::e_inside
                                       : volume_containment::e_outside};
    const point3 p{to_global(l[0], l[1], l[2])};
    EXPECT_EQ(vol.is_inside(p), exp_vc) << "at point " << to_string(p);
  }

  // Inconsistent boundary: The first portal is larger than its face
  const auto large_pt_det = build_rotated_box(host_mr, center, axes, h, 1.5f);
  EXPECT_THROW(detail::check_consistency(large_pt_det), std::invalid_argument);

  // Invalid cuboid: All portals are parallel
  const auto parallel_det =
      build_rotated_box(host_mr, center, axes, h, 1.f, true);
  EXPECT_THROW(detail::check_consistency(parallel_det), std::invalid_argument);
}
