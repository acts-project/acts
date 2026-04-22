// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/type_registry.hpp"

#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <iostream>

namespace detray {

template <bool do_debug = !intersection::contains_pos>
struct select_ray_intersector {
  template <typename mask_t>
  using type = ray_intersector<typename mask_t::shape,
                               typename mask_t::algebra_type, do_debug>;
};

}  // namespace detray

// Test type list implementation
GTEST_TEST(detray_utils, mapped_type_registry) {
  using namespace detray;

  using metadata_t = test::toy_metadata;
  using detector_t = detector<metadata_t>;
  using test_algebra_t = detector_t::algebra_type;
  using mask_types = typename detector_t::masks;
  using mask_id = detector_t::masks::id;

  using mapped_registry_t = types::mapped_registry<
      detector_t::masks, select_ray_intersector<intersection::contains_pos>>;

  types::print<detector_t::masks::type_list>();
  std::cout << std::endl;

  types::print<mapped_registry_t::type_list>();
  std::cout << std::endl;

  //
  // Test the mapping
  //
  constexpr auto idx_array = mapped_registry_t::index_map();

  static_assert(mapped_registry_t::n_types == 3u);
  static_assert(types::size<mapped_registry_t> == 3u);
  static_assert(idx_array.size() == 4u);

  EXPECT_EQ(idx_array[0u], 0u);
  EXPECT_EQ(idx_array[1u], 0u);
  EXPECT_EQ(idx_array[2u], 1u);
  EXPECT_EQ(idx_array[3u], 2u);

  using rectangle_t = types::get<mask_types, mask_id::e_rectangle2D>;
  using trapezoid_t = types::get<mask_types, mask_id::e_trapezoid2D>;
  using cylinder_t = types::get<mask_types, mask_id::e_concentric_cylinder2D>;
  using disc_t = types::get<mask_types, mask_id::e_ring2D>;

  static_assert(types::position<mapped_registry_t, rectangle_t> == 0u);
  static_assert(types::position<mapped_registry_t, trapezoid_t> == 0u);
  static_assert(types::position<mapped_registry_t, cylinder_t> == 1u);
  static_assert(types::position<mapped_registry_t, disc_t> == 2u);

  using cyl_intersector_t =
      types::get<mapped_registry_t, mask_id::e_concentric_cylinder2D>;
  using rect_intersector_t =
      types::get<mapped_registry_t, mask_id::e_rectangle2D>;
  using ring_intersector_t = types::get<mapped_registry_t, mask_id::e_ring2D>;
  using trpz_intersector_t =
      types::get<mapped_registry_t, mask_id::e_trapezoid2D>;

  types::print<types::list<cyl_intersector_t>>();
  std::cout << std::endl;
  types::print<types::list<rect_intersector_t>>();
  std::cout << std::endl;
  types::print<types::list<ring_intersector_t>>();
  std::cout << std::endl;
  types::print<types::list<trpz_intersector_t>>();

  static_assert(
      std::same_as<cyl_intersector_t,
                   ray_intersector<concentric_cylinder2D, test_algebra_t,
                                   intersection::contains_pos>>,
      "Retrieved incorrect type after mapping (cylinder2D)");

  static_assert(std::same_as<rect_intersector_t,
                             ray_intersector<rectangle2D, test_algebra_t,
                                             intersection::contains_pos>>,
                "Retrieved incorrect type after mapping (rectangle2D)");

  static_assert(
      std::same_as<
          ring_intersector_t,
          ray_intersector<ring2D, test_algebra_t, intersection::contains_pos>>,
      "Retrieved incorrect type after mapping (ring2D)");

  static_assert(std::same_as<trpz_intersector_t,
                             ray_intersector<trapezoid2D, test_algebra_t,
                                             intersection::contains_pos>>,
                "Retrieved incorrect type after mapping (trapezoid2D)");

  // Check the indices of the mapped types
  static_assert(types::position<mapped_registry_t, rect_intersector_t> == 0u);
  static_assert(types::position<mapped_registry_t, trpz_intersector_t> == 0u);
  static_assert(types::position<mapped_registry_t, cyl_intersector_t> == 1u);
  static_assert(types::position<mapped_registry_t, ring_intersector_t> == 2u);

  //
  // Test the registry
  //

  // Get ID from the original types
  static_assert(types::id<mapped_registry_t, cylinder_t> ==
                    mask_id::e_concentric_cylinder2D,
                "ID for type cylinder intersector incorrect");
  static_assert(
      types::id<mapped_registry_t, rectangle_t> == mask_id::e_rectangle2D,
      "ID for type rectangle intersector incorrect");
  static_assert(types::id<mapped_registry_t, disc_t> == mask_id::e_ring2D,
                "ID for type ring intersector incorrect");
  static_assert(
      types::id<mapped_registry_t, trapezoid_t> == mask_id::e_trapezoid2D,
      "ID for type trapezoid intersector incorrect");

  static_assert(types::id<mapped_registry_t, cyl_intersector_t> ==
                    mask_id::e_concentric_cylinder2D,
                "ID for type cylinder intersector incorrect");
  static_assert(types::id<mapped_registry_t, rect_intersector_t> ==
                    mask_id::e_rectangle2D,
                "ID for type rectangle intersector incorrect");
  static_assert(
      types::id<mapped_registry_t, ring_intersector_t> == mask_id::e_ring2D,
      "ID for type ring intersector incorrect");
  // From the point of view of the mapped registry, this is the same as for
  // the rectangle shape
  static_assert(types::id<mapped_registry_t, trpz_intersector_t> ==
                    mask_id::e_rectangle2D,
                "ID for type trapezoid intersector incorrect");

  // contains
  static_assert(types::contains<mapped_registry_t, cyl_intersector_t>,
                "'contains' failed for cylinder intersector type");
  static_assert(types::contains<mapped_registry_t, rect_intersector_t>,
                "'contains' failed for rectangle intersector type");
  static_assert(types::contains<mapped_registry_t, ring_intersector_t>,
                "'contains' failed for ring intersector type");
  static_assert(types::contains<mapped_registry_t, trpz_intersector_t>,
                "'contains' failed for trapezoid intersector type");

  static_assert(!types::contains<mapped_registry_t, char>,
                "'contains' failed for 'char' type");
  static_assert(!types::contains<mapped_registry_t, void>,
                "'contains' failed for 'void' type");

  // Is valid
  static_assert(mapped_registry_t::is_valid(mask_id::e_concentric_cylinder2D),
                "ID for cylinder intersector invalid");
  static_assert(mapped_registry_t::is_valid(mask_id::e_rectangle2D),
                "ID for rectangle intersector invalid");
  static_assert(mapped_registry_t::is_valid(mask_id::e_ring2D),
                "ID for ring intersector invalid");
  static_assert(mapped_registry_t::is_valid(mask_id::e_trapezoid2D),
                "ID for trapezoid intersector invalid");
  assert(!mapped_registry_t::is_valid(5u) && "ID '5' not invalid");

  // Get type
  static_assert(
      std::same_as<
          cyl_intersector_t,
          types::get<mapped_registry_t, mask_id::e_concentric_cylinder2D>>,
      "Got incorrect type for 'e_concentric_cylinder2D'");
  static_assert(
      std::same_as<rect_intersector_t,
                   types::get<mapped_registry_t, mask_id::e_rectangle2D>>,
      "Got incorrect type for 'e_rectangle2D'");
  static_assert(std::same_as<ring_intersector_t,
                             types::get<mapped_registry_t, mask_id::e_ring2D>>,
                "Got incorrect type for 'e_ring2D'");
  static_assert(
      std::same_as<trpz_intersector_t,
                   types::get<mapped_registry_t, mask_id::e_trapezoid2D>>,
      "Got incorrect type for 'e_trapezoid2D'");
}
