// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/builders/volume_brute_force_builder.hpp"

#include "detray/builders/detail/bounding_box_creator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/nested_volumes_metadata.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/geometry/volume_descriptor.hpp"
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/common/build_nested_detector.hpp"
#include "detray/test/common/build_telescope_detector.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace detray;

namespace {

using algebra_t = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

/// Collect the daughter volume indices of a volume
struct daughter_collector {
  void operator()(const dindex daughter_idx,
                  std::vector<dindex> &daughters) const {
    daughters.push_back(daughter_idx);
  }
};

/// @returns the sorted daughter indices of volume @p vol_idx
template <typename detector_t>
std::vector<dindex> get_daughters(const detector_t &det, const dindex vol_idx) {
  std::vector<dindex> daughters{};
  tracking_volume{det, vol_idx}
      .template visit_daughter_volumes<daughter_collector>(daughters);
  std::ranges::sort(daughters);

  return daughters;
}

/// Build a detector with two cuboid volumes: An outer volume (index 0) and an
/// inner volume (index 1)
///
/// @param mr memory resource for the detector containers
/// @param outer_mother mother of the outer volume
/// @param inner_mother mother of the inner volume
/// @param inner_center center of the inner volume
///
/// @note The daughter volumes get the default volume brute force accelerator
template <typename metadata_t>
auto build_two_boxes(vecmem::memory_resource &mr, const dindex outer_mother,
                     const dindex inner_mother,
                     const point3 &inner_center = {0.f, 0.f, 0.f}) {
  using builder_t = detector_builder<metadata_t>;
  using detector_t = typename builder_t::detector_type;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  constexpr nav_link_t inv{detail::invalid_value<nav_link_t>()};

  builder_t det_builder{};

  auto *outer = det_builder.new_volume(volume_id::e_cuboid);
  outer->set_mother(outer_mother);
  detail::add_box_portals(*outer, point3{0.f, 0.f, 0.f}, {100.f, 100.f, 100.f},
                          {inv, inv, inv, inv, inv, inv});

  auto *inner = det_builder.new_volume(volume_id::e_cuboid);
  inner->set_mother(inner_mother);
  inner->add_volume_placement(inner_center);
  detail::add_box_portals(*inner, inner_center, {10.f, 10.f, 10.f},
                          {0u, 0u, 0u, 0u, 0u, 0u});

  return det_builder.build(mr, volume_builder_options{});
}

}  // anonymous namespace

/// Unittest: Mother link in the volume descriptor
GTEST_TEST(detray_builders, volume_descriptor_mother) {
  using metadata_t = nested_volumes_metadata<algebra_t>;
  using detector_t = host::detector<metadata_t>;
  using volume_t = typename detector_t::volume_type;

  // The descriptor must stay trivially copyable for device code
  static_assert(std::is_trivially_copyable_v<volume_t>);

  volume_t vol{volume_id::e_cuboid};
  EXPECT_FALSE(vol.has_mother());
  EXPECT_EQ(vol.mother(), dindex_invalid);

  vol.set_mother(3u);
  EXPECT_TRUE(vol.has_mother());
  EXPECT_EQ(vol.mother(), 3u);

  std::stringstream os{};
  os << vol;
  EXPECT_NE(os.str().find("mother = 3"), std::string::npos);

  vol.set_mother(dindex_invalid);
  EXPECT_FALSE(vol.has_mother());
}

/// Unittest: Build the volume hierarchy of the nested test detector
GTEST_TEST(detray_builders, volume_hierarchy_nested_detector) {
  vecmem::host_memory_resource host_mr;

  // Includes the consistency check
  const auto [det, names] = build_nested_detector<algebra_t>(host_mr);

  ASSERT_EQ(det.volumes().size(), 6u);

  // Mother links
  const std::vector<dindex> exp_mothers{dindex_invalid, 0u, 1u, 0u, 0u, 0u};
  for (const auto &vol_desc : det.volumes()) {
    EXPECT_EQ(vol_desc.mother(), exp_mothers.at(vol_desc.index()))
        << names.at(vol_desc.index());
  }

  // Daughter links
  const auto world = tracking_volume{det, 0u};
  EXPECT_TRUE(world.is_top_level());
  EXPECT_TRUE(world.has_daughters());
  EXPECT_EQ(get_daughters(det, 0u), (std::vector<dindex>{1u, 3u, 4u, 5u}));

  const auto box_a = tracking_volume{det, 1u};
  EXPECT_FALSE(box_a.is_top_level());
  EXPECT_EQ(box_a.mother(), 0u);
  EXPECT_TRUE(box_a.has_daughters());
  EXPECT_EQ(get_daughters(det, 1u), (std::vector<dindex>{2u}));

  for (const dindex leaf : {2u, 3u, 4u, 5u}) {
    EXPECT_FALSE(tracking_volume(det, leaf).has_daughters()) << names.at(leaf);
    EXPECT_TRUE(get_daughters(det, leaf).empty()) << names.at(leaf);
  }

  // One daughter list per mother volume (the default volume accelerator is
  // the volume grid, so no index is reserved for the world volume)
  using accel_id = typename std::remove_cvref_t<decltype(det)>::accel::id;
  EXPECT_EQ(
      det.accelerator_store().template size<accel_id::e_volume_brute_force>(),
      2u);
}

/// Unittest: Bounding box of a volume from its portals
GTEST_TEST(detray_builders, volume_hierarchy_bounding_box) {
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = build_nested_detector<algebra_t>(host_mr);

  constexpr scalar tol{1e-3f};

  // World
  auto box = detail::volume_bounding_box(det, 0u);
  auto min_p = box.template loc_min<point3>();
  auto max_p = box.template loc_max<point3>();
  for (unsigned int i = 0u; i < 3u; ++i) {
    EXPECT_NEAR(min_p[i], -200.f, tol);
    EXPECT_NEAR(max_p[i], 200.f, tol);
  }

  // Box C
  box = detail::volume_bounding_box(det, 4u);
  min_p = box.template loc_min<point3>();
  max_p = box.template loc_max<point3>();
  EXPECT_NEAR(min_p[0], -40.f, tol);
  EXPECT_NEAR(max_p[0], 40.f, tol);
  EXPECT_NEAR(min_p[1], 90.f, tol);
  EXPECT_NEAR(max_p[1], 170.f, tol);
  EXPECT_NEAR(min_p[2], -30.f, tol);
  EXPECT_NEAR(max_p[2], 30.f, tol);
}

/// Unittest: Flat detectors do not get a volume hierarchy
GTEST_TEST(detray_builders, volume_hierarchy_flat_detector) {
  vecmem::host_memory_resource host_mr;

  tel_det_config<algebra_t> tel_cfg{
      -20.f * unit<scalar>::mm, -20.f * unit<scalar>::mm,
      20.f * unit<scalar>::mm, 20.f * unit<scalar>::mm};
  tel_cfg.n_surfaces(5u);
  const auto [det, names] =
      build_telescope_detector<algebra_t>(host_mr, tel_cfg);

  for (const auto &vol_desc : det.volumes()) {
    const auto vol = tracking_volume{det, vol_desc};
    EXPECT_TRUE(vol.is_top_level());
    EXPECT_FALSE(vol.has_daughters());
  }

  // Nothing was added to the volume brute force accelerator
  using accel_id = typename std::remove_cvref_t<decltype(det)>::accel::id;
  EXPECT_EQ(
      det.accelerator_store().template size<accel_id::e_volume_brute_force>(),
      0u);
}

/// Unittest: Reserve the world volume lookup, if the brute force
/// accelerator is the default volume accelerator (telescope metadata)
GTEST_TEST(detray_builders, volume_hierarchy_world_lookup) {
  using metadata_t = test::default_telescope_metadata;

  vecmem::host_memory_resource host_mr;
  const auto det = build_two_boxes<metadata_t>(host_mr, dindex_invalid, 0u);

  using accel_id = typename std::remove_cvref_t<decltype(det)>::accel::id;
  static_assert(accel_id::e_volume_brute_force == accel_id::e_volume_default);

  // Index 0: top-level volumes, index 1: daughters of volume 0
  const auto &bf_coll =
      det.accelerator_store().template get<accel_id::e_volume_brute_force>();
  ASSERT_EQ(bf_coll.size(), 2u);

  std::vector<dindex> top_level{};
  for (const dindex idx : bf_coll[0].all()) {
    top_level.push_back(idx);
  }
  EXPECT_EQ(top_level, (std::vector<dindex>{0u}));

  EXPECT_EQ(det.volume(0u)
                .template accel_link<
                    std::remove_cvref_t<decltype(det)>::geo_obj_ids::e_volume>()
                .index(),
            1u);
  EXPECT_EQ(get_daughters(det, 0u), (std::vector<dindex>{1u}));

  // Point lookup through the world volume list and the daughter volumes
  EXPECT_EQ(det.volume_index(point3{1.f, 2.f, 3.f}), 1u);
  EXPECT_EQ(det.volume_index(point3{50.f, 2.f, 3.f}), 0u);
  EXPECT_EQ(det.volume_index(point3{150.f, 2.f, 3.f}), dindex_invalid);
  // Outside of the detector: Invalid volume descriptor
  EXPECT_EQ(det.volume(point3{150.f, 2.f, 3.f}).index(), dindex_invalid);
  EXPECT_EQ(det.volume(point3{50.f, 2.f, 3.f}).index(), 0u);
}

/// Unittest: The world root can only be registered with one builder, which
/// must be the first volume brute force builder
GTEST_TEST(detray_builders, volume_brute_force_world_root) {
  using metadata_t = test::default_telescope_metadata;
  using builder_t = detector_builder<metadata_t>;
  using detector_t = typename builder_t::detector_type;
  using bf_builder_t = volume_brute_force_builder<detector_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  constexpr nav_link_t inv{detail::invalid_value<nav_link_t>()};

  vecmem::host_memory_resource host_mr;

  // Outer box (top-level) with an inner box
  auto add_boxes = [](builder_t &det_builder) {
    auto *outer = det_builder.new_volume(volume_id::e_cuboid);
    detail::add_box_portals(*outer, point3{0.f, 0.f, 0.f},
                            {100.f, 100.f, 100.f},
                            {inv, inv, inv, inv, inv, inv});
    auto *inner = det_builder.new_volume(volume_id::e_cuboid);
    inner->set_mother(0u);
    detail::add_box_portals(*inner, point3{0.f, 0.f, 0.f}, {10.f, 10.f, 10.f},
                            {0u, 0u, 0u, 0u, 0u, 0u});
  };

  // World root registered twice
  builder_t twice_builder{};
  add_boxes(twice_builder);
  twice_builder.add_volume_accelerator<bf_builder_t>()->add_root(dindex_invalid,
                                                                 {0u});
  twice_builder.add_volume_accelerator<bf_builder_t>()->add_root(dindex_invalid,
                                                                 {0u});
  EXPECT_THROW(twice_builder.build(host_mr, volume_builder_options{}),
               std::invalid_argument);

  // The first builder already built the world volume lookup
  builder_t late_builder{};
  add_boxes(late_builder);
  late_builder.add_volume_accelerator<bf_builder_t>()->add_root(0u, {1u});
  late_builder.add_volume_accelerator<bf_builder_t>()->add_root(dindex_invalid,
                                                                {0u});
  EXPECT_THROW(late_builder.build(host_mr, volume_builder_options{}),
               std::invalid_argument);

  // Valid: Explicit world root and daughters in one builder
  builder_t det_builder{};
  add_boxes(det_builder);
  auto *bf_builder = det_builder.add_volume_accelerator<bf_builder_t>();
  bf_builder->add_root(dindex_invalid, {0u});
  bf_builder->add_root(0u, {1u});
  const auto det = det_builder.build(host_mr, volume_builder_options{});
  EXPECT_NO_THROW(detail::check_consistency(det));
  EXPECT_EQ(det.volume_index(point3{1.f, 2.f, 3.f}), 1u);
}

/// Unittest: Broken volume hierarchies are rejected by the builder (cannot
/// build) or by the consistency checker (the navigation would fail)
GTEST_TEST(detray_builders, volume_hierarchy_errors) {
  using metadata_t = nested_volumes_metadata<algebra_t>;

  vecmem::host_memory_resource host_mr;

  auto check = [](const auto &det) { detail::check_consistency(det); };

  // Valid hierarchy
  EXPECT_NO_THROW(
      check(build_two_boxes<metadata_t>(host_mr, dindex_invalid, 0u)));

  // Build errors: Mother index out of bounds (no root volume to link)
  EXPECT_THROW(build_two_boxes<metadata_t>(host_mr, dindex_invalid, 2u),
               std::invalid_argument);

  // The detector type cannot hold daughter volumes: The daughter is not part
  // of any volume accelerator (only a warning, until volume grids exist)
  EXPECT_NO_THROW(
      check(build_two_boxes<test::toy_metadata>(host_mr, dindex_invalid, 0u)));

  // Consistency errors: Mother index of an ancestor out of bounds (the toy
  // detector type has no volume accelerator, so the builder does not fail)
  EXPECT_THROW(check(build_two_boxes<test::toy_metadata>(host_mr, 1u, 500u)),
               std::invalid_argument);

  // Consistency errors: Volume is its own mother
  EXPECT_THROW(check(build_two_boxes<metadata_t>(host_mr, dindex_invalid, 1u)),
               std::invalid_argument);

  // Consistency errors: Cycle (the volumes are mothers of each other)
  EXPECT_THROW(check(build_two_boxes<metadata_t>(host_mr, 1u, 0u)),
               std::invalid_argument);

  // Consistency errors: Daughter is not contained in the mother
  const auto outside_det = build_two_boxes<metadata_t>(
      host_mr, dindex_invalid, 0u, point3{95.f, 0.f, 0.f});
  EXPECT_THROW(check(outside_det), std::invalid_argument);
}

/// Unittest: Root volumes are registered explicitly with a volume accelerator
/// builder, all other mother volumes get the default brute force builder
GTEST_TEST(detray_builders, volume_brute_force_builder) {
  using metadata_t = nested_volumes_metadata<algebra_t>;
  using builder_t = detector_builder<metadata_t>;
  using detector_t = typename builder_t::detector_type;
  using bf_builder_t = volume_brute_force_builder<detector_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  constexpr nav_link_t inv{detail::invalid_value<nav_link_t>()};

  vecmem::host_memory_resource host_mr;

  // Box 0 contains the boxes 1 and 2, box 1 contains box 3
  auto add_boxes = [](builder_t &det_builder) {
    auto add_box = [&det_builder](const dindex mother, const point3 &center,
                                  const scalar half_length,
                                  const nav_link_t link) {
      auto *v_builder = det_builder.new_volume(volume_id::e_cuboid);
      v_builder->set_mother(mother);
      v_builder->add_volume_placement(center);
      detail::add_box_portals(*v_builder, center,
                              {half_length, half_length, half_length},
                              {link, link, link, link, link, link});
    };
    add_box(dindex_invalid, {0.f, 0.f, 0.f}, 100.f, inv);
    add_box(0u, {-50.f, 0.f, 0.f}, 40.f, 0u);
    add_box(0u, {50.f, 0.f, 0.f}, 40.f, 0u);
    add_box(1u, {-50.f, 0.f, 0.f}, 10.f, 1u);
  };

  // Builder specific checks
  bf_builder_t bf_builder{};
  EXPECT_EQ(bf_builder.name(), "volume_brute_force_builder");
  EXPECT_THROW(bf_builder.add_root(0u, {}), std::invalid_argument);
  EXPECT_THROW(bf_builder.add_root(0u, {1u, 1u}), std::invalid_argument);
  // The brute force accelerator is not the default volume accelerator
  EXPECT_THROW(bf_builder.add_root(dindex_invalid, {0u}),
               std::invalid_argument);
  bf_builder.add_root(0u, {1u, 2u});
  EXPECT_THROW(bf_builder.add_root(0u, {1u}), std::invalid_argument);
  EXPECT_EQ(bf_builder.roots(), (std::vector<dindex>{0u}));

  // Explicit root: Box 0; Default builder: Box 1
  builder_t det_builder{};
  add_boxes(det_builder);
  auto *v_acc_builder = det_builder.add_volume_accelerator<bf_builder_t>();
  ASSERT_NE(v_acc_builder, nullptr);
  v_acc_builder->add_root(0u, {2u, 1u});
  EXPECT_EQ(det_builder.n_volume_accelerators(), 1u);

  const auto det = det_builder.build(host_mr, volume_builder_options{});
  EXPECT_NO_THROW(detail::check_consistency(det));

  EXPECT_EQ(get_daughters(det, 0u), (std::vector<dindex>{1u, 2u}));
  EXPECT_EQ(get_daughters(det, 1u), (std::vector<dindex>{3u}));

  // The explicit builder runs before the default builder
  constexpr auto vol_obj_id{detector_t::geo_obj_ids::e_volume};
  EXPECT_EQ(det.volume(0u).template accel_link<vol_obj_id>().index(), 0u);
  EXPECT_EQ(det.volume(1u).template accel_link<vol_obj_id>().index(), 1u);

  // A daughter that is missing in the explicit list is not part of any
  // volume accelerator (only a warning, until volume grids exist)
  builder_t missing_builder{};
  add_boxes(missing_builder);
  missing_builder.add_volume_accelerator<bf_builder_t>()->add_root(0u, {1u});
  const auto missing_det =
      missing_builder.build(host_mr, volume_builder_options{});
  EXPECT_NO_THROW(detail::check_consistency(missing_det));

  // The explicit daughter list must match the mother links
  builder_t wrong_builder{};
  add_boxes(wrong_builder);
  wrong_builder.add_volume_accelerator<bf_builder_t>()->add_root(0u, {1u, 2u});
  wrong_builder.add_volume_accelerator<bf_builder_t>()->add_root(1u, {3u, 2u});
  const auto wrong_det = wrong_builder.build(host_mr, volume_builder_options{});
  EXPECT_THROW(detail::check_consistency(wrong_det), std::invalid_argument);

  // A root must not be registered with two builders
  builder_t twice_builder{};
  add_boxes(twice_builder);
  twice_builder.add_volume_accelerator<bf_builder_t>()->add_root(0u, {1u, 2u});
  twice_builder.add_volume_accelerator<bf_builder_t>()->add_root(0u, {1u, 2u});
  EXPECT_THROW(twice_builder.build(host_mr, volume_builder_options{}),
               std::invalid_argument);
}
