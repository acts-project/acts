// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/geometry/surface.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/utils/geometry_utils.hpp"  //< cos incidence angle

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

namespace {

/// Define mask types
enum class mask_id : unsigned int {
  e_unmasked = 0u,
};

std::ostream& operator<<(std::ostream& os, mask_id /*mid*/) {
  os << "e_unmasked";
  return os;
}

/// Define material types
enum class material_id : unsigned int {
  e_material_slab = 0u,
};

std::ostream& operator<<(std::ostream& os, material_id /*mid*/) {
  os << "e_material_slab";
  return os;
}

constexpr detray::test::scalar tol{5e-5f};

detray::dvector3D<detray::test::algebra> test_dir{0.f, 0.f, 1.f};

}  // anonymous namespace

// This tests the construction of a surface descriptor object
GTEST_TEST(detray_geometry, surface_descriptor) {
  using namespace detray;

  using mask_link_t = dtyped_index<mask_id, dindex>;
  using material_link_t = dtyped_index<material_id, dindex>;

  mask_link_t mask_id{mask_id::e_unmasked, 0u};
  material_link_t material_id{material_id::e_material_slab, 0u};

  surface_descriptor<mask_link_t, material_link_t> desc(
      1u, mask_id, material_id, 2u, surface_id::e_sensitive);

  static_assert(sizeof(decltype(desc)) == 16);

  // Test access
  ASSERT_EQ(desc.transform(), 1u);
  ASSERT_EQ(desc.volume(), 2u);
  ASSERT_EQ(desc.id(), surface_id::e_sensitive);
  ASSERT_FALSE(desc.is_portal());
  ASSERT_FALSE(desc.is_passive());
  ASSERT_TRUE(desc.is_sensitive());
  ASSERT_EQ(desc.mask(), mask_id);
  ASSERT_EQ(desc.material(), material_id);

  // Test setters
  desc.set_volume(5u);
  desc.set_id(surface_id::e_portal);
  desc.set_index(6u);
  desc.update_transform(7u);
  desc.update_mask(7u);
  desc.update_material(7u);

  ASSERT_EQ(desc.transform(), 8u);
  ASSERT_EQ(desc.volume(), 5u);
  ASSERT_EQ(desc.id(), surface_id::e_portal);
  ASSERT_EQ(desc.index(), 6u);
  ASSERT_TRUE(desc.is_portal());
  ASSERT_FALSE(desc.is_passive());
  ASSERT_FALSE(desc.is_sensitive());
  ASSERT_EQ(desc.mask().id(), mask_id::e_unmasked);
  ASSERT_EQ(desc.mask().index(), 7u);
  ASSERT_EQ(desc.material().id(), material_id::e_material_slab);
  ASSERT_EQ(desc.material().index(), 7u);
}

/// This tests the functionality of a tracking surface interface for the toy
/// detector
GTEST_TEST(detray_geometry, surface_toy_detector) {
  using namespace detray;

  using detector_t = detector<test::toy_metadata>;

  using test_algebra = detector_t::algebra_type;
  using scalar = geometry::surface<detector_t>::scalar_type;
  using point2 = geometry::surface<detector_t>::point2_type;
  using point3 = geometry::surface<detector_t>::point3_type;
  using vector3 = geometry::surface<detector_t>::vector3_type;

  vecmem::host_memory_resource host_mr;
  toy_det_config<test::scalar> toy_cfg{};
  toy_cfg.use_material_maps(true).do_check(true);
  const auto [toy_det, names] =
      build_toy_detector<test_algebra>(host_mr, toy_cfg);

  auto ctx = typename detector_t::geometry_context{};

  //
  // Disc
  //
  const auto disc_descr = toy_det.surfaces()[1u];
  const auto disc = geometry::surface{toy_det, disc_descr};

  // IDs
  ASSERT_EQ(disc.identifier(), disc_descr.identifier());
  ASSERT_EQ(disc.volume(), 0u);
  ASSERT_EQ(disc.index(), 1u);
  ASSERT_EQ(disc.id(), surface_id::e_portal);
  ASSERT_EQ(disc.shape_id(), detector_t::masks::id::e_ring2D);
  ASSERT_EQ(disc.shape_name(), "ring2D");
  ASSERT_FALSE(disc.is_sensitive());
  ASSERT_FALSE(disc.is_passive());
  ASSERT_TRUE(disc.is_portal());

  // Transformation matrix
  const auto disc_translation =
      disc.transform(ctx).translation();  // beampipe portal
  ASSERT_NEAR(disc_translation[0], 0.f, tol);
  ASSERT_NEAR(disc_translation[1], 0.f, tol);
  ASSERT_NEAR(disc_translation[2], -827.5f, tol);
  auto center = disc.center(ctx);
  ASSERT_NEAR(center[0], 0.f, tol);
  ASSERT_NEAR(center[1], 0.f, tol);
  ASSERT_NEAR(center[2], -827.5f, tol);

  // Surface normal
  const auto z_axis = vector3{0.f, 0.f, 1.f};
  const vector3 normal_3D = disc.normal(ctx, point3{0.f, 0.f, 0.f});
  const vector3 normal_2D = disc.normal(ctx, point2{0.f, 0.f});
  // trigger all code paths
  ASSERT_NEAR(normal_3D[0], z_axis[0], tol);
  ASSERT_NEAR(normal_3D[1], z_axis[1], tol);
  ASSERT_NEAR(normal_3D[2], z_axis[2], tol);
  ASSERT_NEAR(normal_2D[0], z_axis[0], tol);
  ASSERT_NEAR(normal_2D[1], z_axis[1], tol);
  ASSERT_NEAR(normal_2D[2], z_axis[2], tol);

  // Cos incidence angle
  vector3 dir = vector::normalize(vector3{1.f, 0.f, 1.f});
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point3{0.f, 0.f, 0.f}),
              constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point2{0.f, 0.f}),
              constant<scalar>::inv_sqrt2, tol);

  dir = vector::normalize(vector3{1.f, 1.f, 0.f});
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point3{1.f, 0.5f, 0.f}), 0.f, tol);
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point2{1.f, 0.5f}), 0.f, tol);

  dir = vector::normalize(vector3{static_cast<scalar>(1.f),
                                  static_cast<scalar>(0.f),
                                  constant<scalar>::pi});
  scalar cos_inc_angle{constant<scalar>::pi /
                       std::sqrt(1.f + std::pow(constant<scalar>::pi, 2.f))};
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point3{2.f, 1.f, 0.f}), cos_inc_angle,
              tol);
  ASSERT_NEAR(cos_angle(ctx, disc, dir, point2{2.f, 1.f}), cos_inc_angle, tol);

  // Coordinate transformations
  point3 glob_pos = {4.f, 7.f, 4.f};
  point3 local = disc.global_to_local(ctx, glob_pos, test_dir);
  point2 bound = disc.global_to_bound(ctx, glob_pos, test_dir);

  ASSERT_NEAR(local[0], std::sqrt(65.f), tol);
  ASSERT_NEAR(local[1], std::atan2(7.f, 4.f), tol);
  ASSERT_NEAR(bound[0], local[0], tol);
  ASSERT_NEAR(bound[1], local[1], tol);

  // Roundtrip
  point3 global = disc.local_to_global(ctx, local, test_dir);
  point3 global2 = disc.local_to_global(ctx, bound, test_dir);

  ASSERT_NEAR(glob_pos[0], global[0], tol);
  ASSERT_NEAR(glob_pos[1], global[1], tol);
  ASSERT_NEAR(glob_pos[2], global[2], tol);

  ASSERT_NEAR(global2[0], glob_pos[0], tol);
  ASSERT_NEAR(global2[1], glob_pos[1], tol);
  // The bound transform assumes the point is on surface
  ASSERT_NEAR(global2[2], disc_translation[2], tol);

  // Test the material
  ASSERT_TRUE(disc.has_material());
  const auto* mat_param = disc.material_parameters({0.f, 0.f});
  ASSERT_TRUE(mat_param);
  ASSERT_EQ(*mat_param, toy_cfg.mapped_material());

  //
  // Rectangle
  //
  const auto rec_descr = toy_det.surfaces()[586u];
  const auto rec = geometry::surface{toy_det, rec_descr};

  // IDs
  ASSERT_EQ(rec.identifier(), rec_descr.identifier());
  ASSERT_EQ(rec.volume(), 9u);
  ASSERT_EQ(rec.index(), 586u);
  ASSERT_EQ(rec.id(), surface_id::e_sensitive);
  ASSERT_EQ(rec.shape_id(), detector_t::masks::id::e_rectangle2D);
  ASSERT_TRUE(rec.is_sensitive());
  ASSERT_FALSE(rec.is_passive());
  ASSERT_FALSE(rec.is_portal());

  // Transformation matrix
  const auto& rec_translation = rec.transform(ctx).translation();
  ASSERT_NEAR(rec_translation[0], -71.902099f, tol);
  ASSERT_NEAR(rec_translation[1], -7.081735f, tol);
  ASSERT_NEAR(rec_translation[2], -455.f, tol);
  center = rec.center(ctx);
  ASSERT_NEAR(center[0], -71.902099f, tol);
  ASSERT_NEAR(center[1], -7.081735f, tol);
  ASSERT_NEAR(center[2], -455.f, tol);

  // Surface normal
  // trigger all code paths
  global = rec.transform(ctx).vector_to_global(z_axis);
  ASSERT_EQ(rec.normal(ctx, point3{0.f, 0.f, 0.f}), global);
  ASSERT_EQ(rec.normal(ctx, point2{0.f, 0.f}), global);

  // Incidence angle
  dir = vector::normalize(global);
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point3{0.f, 0.f, 0.f}), 1.f, tol);
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point2{0.f, 0.f}), 1.f, tol);

  dir = vector::normalize(
      vector3{static_cast<scalar>(0.f), -global[2], global[1]});
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point3{1.f, 0.5f, 0.f}), 0.f, tol);
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point2{1.f, 0.5f}), 0.f, tol);

  dir = vector::normalize(vector3{static_cast<scalar>(1.f),
                                  static_cast<scalar>(0.f),
                                  constant<scalar>::pi});
  cos_inc_angle = std::fabs(vector::dot(dir, global));
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point3{2.f, 1.f, 0.f}), cos_inc_angle,
              tol);
  ASSERT_NEAR(cos_angle(ctx, rec, dir, point2{2.f, 1.f}), cos_inc_angle, tol);

  // Coordinate transformation roundtrip
  glob_pos = {4.f, 7.f, 4.f};

  local = rec.global_to_local(ctx, glob_pos, test_dir);
  global = rec.local_to_global(ctx, local, test_dir);
  ASSERT_NEAR(glob_pos[0], global[0], tol);
  ASSERT_NEAR(glob_pos[1], global[1], tol);
  ASSERT_NEAR(glob_pos[2], global[2], tol);

  glob_pos = {-71.902099f, -7.081735f, -460.f};

  local = rec.global_to_local(ctx, glob_pos, test_dir);
  global = rec.local_to_global(ctx, local, test_dir);
  ASSERT_NEAR(glob_pos[0], global[0], tol);
  ASSERT_NEAR(glob_pos[1], global[1], tol);
  ASSERT_NEAR(glob_pos[2], global[2], tol);

  bound = rec.global_to_bound(ctx, glob_pos, test_dir);
  global = rec.local_to_global(ctx, bound, test_dir);
  ASSERT_NEAR(global[0], glob_pos[0], tol);
  ASSERT_NEAR(global[1], glob_pos[1], tol);
  ASSERT_NEAR(global[2], glob_pos[2], tol);

  // Test the material (no material on sensitive surfaces)
  ASSERT_FALSE(rec.has_material());

  //
  // Concentric Cylinder
  //
  const auto cyl_descr = toy_det.surfaces()[3u];
  const auto cyl = geometry::surface{toy_det, cyl_descr};

  // IDs
  ASSERT_EQ(cyl.identifier(), cyl_descr.identifier());
  ASSERT_EQ(cyl.volume(), 0u);
  ASSERT_EQ(cyl.index(), 3u);
  ASSERT_EQ(cyl.id(), surface_id::e_portal);
  ASSERT_EQ(cyl.shape_id(), detector_t::masks::id::e_concentric_cylinder2D);
  ASSERT_FALSE(cyl.is_sensitive());
  ASSERT_FALSE(cyl.is_passive());
  ASSERT_TRUE(cyl.is_portal());

  // Transformation matrix
  const auto& cyl_translation = cyl.transform(ctx).translation();
  ASSERT_NEAR(cyl_translation[0], 0.f, tol);
  ASSERT_NEAR(cyl_translation[1], 0.f, tol);
  ASSERT_NEAR(cyl_translation[2], 0.f, tol);
  center = cyl.center(ctx);
  ASSERT_NEAR(center[0], 0.f, tol);
  ASSERT_NEAR(center[1], 0.f, tol);
  ASSERT_NEAR(center[2], 0.f, tol);

  // Surface normal
  // trigger all code paths
  constexpr scalar r{25.f * unit<scalar>::mm};
  const vector3 x_axis{1.f, 0.f, 0.f};
  ASSERT_NEAR(
      vector::norm(cyl.normal(ctx, point3{static_cast<scalar>(0.f),
                                          static_cast<scalar>(0.f), r}) -
                   x_axis),
      0.f, tol);
  ASSERT_NEAR(vector::norm(cyl.normal(ctx, point2{0.f, 0.f}) - x_axis),
              static_cast<scalar>(0.f), tol);
  ASSERT_NEAR(
      vector::norm(cyl.normal(ctx, point3{r * constant<scalar>::pi,
                                          static_cast<scalar>(0.f), r}) +
                   x_axis),
      0.f, tol);
  ASSERT_NEAR(vector::norm(cyl.normal(ctx, point2{r * constant<scalar>::pi,
                                                  static_cast<scalar>(0.f)}) +
                           x_axis),
              0.f, tol);

  const vector3 y_axis{0.f, 1.f, 0.f};
  ASSERT_NEAR(
      vector::norm(cyl.normal(ctx, point3{r * constant<scalar>::pi_2,
                                          static_cast<scalar>(0.f), r}) -
                   y_axis),
      0.f, tol);
  ASSERT_NEAR(vector::norm(cyl.normal(ctx, point2{r * constant<scalar>::pi_2,
                                                  static_cast<scalar>(0.f)}) -
                           y_axis),
              0.f, tol);

  // Incidence angle
  ASSERT_NEAR(
      cos_angle(ctx, cyl, x_axis,
                point3{r * constant<scalar>::pi, static_cast<scalar>(0.f), r}),
      1.f, tol);
  ASSERT_NEAR(
      cos_angle(ctx, cyl, x_axis,
                point2{r * constant<scalar>::pi, static_cast<scalar>(0.f)}),
      1.f, tol);

  ASSERT_NEAR(cos_angle(ctx, cyl, z_axis,
                        point3{r * constant<scalar>::pi_2,
                               static_cast<scalar>(0.f), r}),
              0.f, tol);
  ASSERT_NEAR(
      cos_angle(ctx, cyl, z_axis,
                point2{r * constant<scalar>::pi_2, static_cast<scalar>(0.f)}),
      0.f, tol);

  dir = vector::normalize(vector3{1.f, 1.f, 0.f});
  ASSERT_NEAR(
      cos_angle(ctx, cyl, dir,
                point3{static_cast<scalar>(0.f), static_cast<scalar>(1.f), r}),
      constant<scalar>::inv_sqrt2, tol);
  ASSERT_NEAR(cos_angle(ctx, cyl, dir, point2{0.f, 1.f}),
              constant<scalar>::inv_sqrt2, tol);

  // Coordinate transformation roundtrip
  glob_pos = {4.f, 7.f, 4.f};

  local = cyl.global_to_local(ctx, glob_pos, test_dir);
  global = cyl.local_to_global(ctx, local, test_dir);
  ASSERT_NEAR(glob_pos[0], global[0], tol);
  ASSERT_NEAR(glob_pos[1], global[1], tol);
  ASSERT_NEAR(glob_pos[2], global[2], tol);

  glob_pos = {constant<scalar>::inv_sqrt2 * r, constant<scalar>::inv_sqrt2 * r,
              static_cast<scalar>(2.f)};

  local = cyl.global_to_local(ctx, glob_pos, test_dir);
  global = cyl.local_to_global(ctx, local, test_dir);
  ASSERT_NEAR(glob_pos[0], global[0], tol);
  ASSERT_NEAR(glob_pos[1], global[1], tol);
  ASSERT_NEAR(glob_pos[2], global[2], tol);

  bound = cyl.global_to_bound(ctx, glob_pos, test_dir);
  global = cyl.local_to_global(ctx, bound, test_dir);
  ASSERT_NEAR(global[0], glob_pos[0], tol);
  ASSERT_NEAR(global[1], glob_pos[1], tol);
  ASSERT_NEAR(global[2], glob_pos[2], tol);

  // Test the material
  ASSERT_TRUE(cyl.has_material());
  mat_param = cyl.material_parameters({0.f, 0.f});
  ASSERT_TRUE(mat_param);
  ASSERT_EQ(*mat_param, toy_cfg.mapped_material());
}

/// This tests the functionality of a tracking surface interface for a wire
/// chamber detector
GTEST_TEST(detray_geometry, surface_wire_chamber) {
  using namespace detray;

  using metadata_t = test::wire_chamber_metadata;
  using detector_t = detector<metadata_t>;

  using test_algebra = metadata_t::algebra_type;
  using scalar = geometry::surface<detector_t>::scalar_type;
  using point2 = geometry::surface<detector_t>::point2_type;
  using point3 = geometry::surface<detector_t>::point3_type;
  using vector3 = geometry::surface<detector_t>::vector3_type;

  vecmem::host_memory_resource host_mr;
  wire_chamber_config<scalar> cfg{};
  const auto [wire_chmbr, names] =
      build_wire_chamber<test_algebra>(host_mr, cfg);

  auto ctx = typename detector_t::geometry_context{};

  //
  // Line
  //
  const auto line_descr = wire_chmbr.surfaces()[23u];
  const auto line = geometry::surface{wire_chmbr, line_descr};

  // IDs
  ASSERT_EQ(line.identifier(), line_descr.identifier());
  ASSERT_EQ(line.volume(), 1u);
  ASSERT_EQ(line.index(), 23u);
  ASSERT_EQ(line.id(), surface_id::e_sensitive);
  ASSERT_EQ(line.shape_id(), detector_t::masks::id::e_drift_cell);
  ASSERT_TRUE(line.is_sensitive());
  ASSERT_FALSE(line.is_passive());
  ASSERT_FALSE(line.is_portal());

  // Transformation matrix
  const auto line_translation = line.transform(ctx).translation();
  ASSERT_NEAR(line_translation[0], 412.858582f, tol);
  ASSERT_NEAR(line_translation[1], 299.412414f, tol);
  ASSERT_NEAR(line_translation[2], 0.f, tol);
  auto center = line.center(ctx);
  ASSERT_NEAR(center[0], 412.858582f, tol);
  ASSERT_NEAR(center[1], 299.412414f, tol);
  ASSERT_NEAR(center[2], 0.f, tol);

  // Surface normal
  const auto z_axis = vector3{0.f, 0.f, 1.f};
  point3 global = line.transform(ctx).vector_to_global(z_axis);
  const vector3 normal_3D = line.normal(ctx, point3{1.f, 0.f, 0.f});
  const vector3 normal_2D = line.normal(ctx, point2{1.f, 0.f});
  const vector3 neg_normal_3D = line.normal(ctx, point3{-1.f, 0.f, 0.f});
  const vector3 neg_normal_2D = line.normal(ctx, point2{-1.f, 0.f});
  // trigger all code paths
  ASSERT_NEAR(normal_3D[0], global[0], tol);
  ASSERT_NEAR(normal_3D[1], global[1], tol);
  ASSERT_NEAR(normal_3D[2], global[2], tol);
  ASSERT_NEAR(normal_2D[0], global[0], tol);
  ASSERT_NEAR(normal_2D[1], global[1], tol);
  ASSERT_NEAR(normal_2D[2], global[2], tol);
  ASSERT_NEAR(neg_normal_3D[0], global[0], tol);
  ASSERT_NEAR(neg_normal_3D[1], global[1], tol);
  ASSERT_NEAR(neg_normal_3D[2], global[2], tol);
  ASSERT_NEAR(neg_normal_2D[0], global[0], tol);
  ASSERT_NEAR(neg_normal_2D[1], global[1], tol);
  ASSERT_NEAR(neg_normal_2D[2], global[2], tol);

  // Cos incidence angle
  vector3 dir = vector::normalize(global);
  ASSERT_NEAR(cos_angle(ctx, line, dir, point3{1.f, 0.f, 0.f}), 1.f, tol);
  ASSERT_NEAR(cos_angle(ctx, line, dir, point2{1.f, 0.f}), 1.f, tol);

  dir = vector::normalize(
      vector3{static_cast<scalar>(0.f), -global[2], global[1]});
  ASSERT_NEAR(cos_angle(ctx, line, dir, point3{1.f, 100.f, 0.f}), 0.f, tol);
  ASSERT_NEAR(cos_angle(ctx, line, dir, point2{1.f, 100.f}), 0.f, tol);

  dir = vector3{-0.685475f, -0.0404595f, 0.726971f};
  ASSERT_NEAR(cos_angle(ctx, line, dir, point3{2.f, 1.f, 0.f}),
              constant<scalar>::inv_sqrt2, 0.0005);
  ASSERT_NEAR(cos_angle(ctx, line, dir, point2{2.f, 1.f}),
              constant<scalar>::inv_sqrt2, 0.0005);

  // Coordinate transformation roundtrip
  point3 glob_pos = {4.f, 7.f, 4.f};

  point3 local = line.global_to_local(ctx, glob_pos, test_dir);
  global = line.local_to_global(ctx, local, dir);

  // @TODO: Needs a reduced tolerance, why?
  scalar red_tol{0.00015f};
  ASSERT_NEAR(glob_pos[0], global[0], red_tol);
  ASSERT_NEAR(glob_pos[1], global[1], red_tol);
  ASSERT_NEAR(glob_pos[2], global[2], red_tol);

  glob_pos = center;

  local = line.global_to_local(ctx, glob_pos, test_dir);
  global = line.local_to_global(ctx, local, test_dir);
  red_tol = 7.f * 1e-5f;
  ASSERT_NEAR(glob_pos[0], global[0], red_tol);
  ASSERT_NEAR(glob_pos[1], global[1], red_tol);
  ASSERT_NEAR(glob_pos[2], global[2], red_tol);

  point2 bound = line.global_to_bound(ctx, glob_pos, test_dir);
  global = line.local_to_global(ctx, bound, test_dir);
  ASSERT_NEAR(global[0], glob_pos[0], red_tol);
  ASSERT_NEAR(global[1], glob_pos[1], red_tol);
  ASSERT_NEAR(global[2], glob_pos[2], red_tol);

  // Test the material
  ASSERT_TRUE(line.has_material());
  const auto* mat_param = line.material_parameters({0.f, 0.f});
  ASSERT_TRUE(mat_param);
  ASSERT_EQ(*mat_param, tungsten<scalar>());
}
