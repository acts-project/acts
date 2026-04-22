// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/indexing.hpp"

#include "detray/utils/concepts.hpp"

// Google test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <iostream>
#include <limits>

using namespace detray;

namespace {

/// Define mask types
enum class mask_id : std::uint8_t {
  e_unmasked = 0u,
};

std::ostream& operator<<(std::ostream& os, mask_id /*mid*/) {
  os << "e_unmasked";
  return os;
}

}  // namespace

/// Test the index range
GTEST_TEST(detray_core, index_range) {
  using index_t = detail::index_range<unsigned int>;
  constexpr auto inv_idx{std::numeric_limits<std::uint32_t>::max()};

  static_assert(concepts::interval<index_t>);

  // Check a empty index range (invalid)
  auto ir = index_t{};
  EXPECT_EQ(ir.lower(), inv_idx);
  EXPECT_EQ(ir.upper(), inv_idx);
  EXPECT_TRUE(ir.is_invalid());

  // Check the values after setting them
  ir.set_lower(5u).set_upper(42u);

  EXPECT_EQ(ir.lower(), 5u);
  EXPECT_EQ(ir.upper(), 42u);
  EXPECT_EQ(detail::get<0>(ir), 5u);
  EXPECT_EQ(detail::get<1>(ir), 42u);
  EXPECT_EQ(ir.size(), 37u);

  // Check invalid link
  EXPECT_FALSE(ir.is_invalid());
  ir.set_lower(inv_idx);
  EXPECT_TRUE(ir.is_invalid());
  ir.set_lower(10u);
  EXPECT_FALSE(ir.is_invalid());

  ir.set_upper(inv_idx);
  EXPECT_TRUE(ir.is_invalid());
  ir.set_upper(20u);
  EXPECT_FALSE(ir.is_invalid());

  EXPECT_EQ(ir.lower(), 10u);
  EXPECT_EQ(ir.upper(), 20u);
  EXPECT_EQ(detail::get<0>(ir), 10u);
  EXPECT_EQ(detail::get<1>(ir), 20u);
  EXPECT_EQ(ir.upper(), 20u);
  EXPECT_EQ(ir.size(), 10u);

  // Test equality operator
  EXPECT_TRUE(ir == ir);
  index_t other = ir;
  EXPECT_TRUE(ir == other);
  ir.set_lower(13u);
  EXPECT_FALSE(ir == other);

  // Test print operator
  std::clog << ir << std::endl;

  // Test arithmetic operators
  index_t ir2 = ir + index_t{3u, 5u};

  EXPECT_EQ(ir2.lower(), 3u);
  EXPECT_EQ(ir2.upper(), 20u);
  EXPECT_EQ(ir2.size(), 17u);
  ir2 = ir2 + 25u;
  EXPECT_EQ(ir2.lower(), 3u);
  EXPECT_EQ(ir2.upper(), 45u);
  EXPECT_EQ(ir2.size(), 42u);

  index_t ir3 = ir - index_t{3u, 15u};

  EXPECT_EQ(ir3.lower(), 13u);
  EXPECT_EQ(ir3.upper(), 15u);
  EXPECT_EQ(ir3.size(), 2u);
  ir3 = ir3 + 2u;
  EXPECT_EQ(ir3.lower(), 13u);
  EXPECT_EQ(ir3.upper(), 17u);
  EXPECT_EQ(ir3.size(), 4u);

  // Comparison operators
  EXPECT_TRUE(ir3 < ir2);
  EXPECT_TRUE(ir3 != ir2);
  EXPECT_TRUE(ir2 > ir3);

  EXPECT_FALSE(ir3 > ir2);
  EXPECT_FALSE(ir3 == ir2);
  EXPECT_FALSE(ir2 < ir3);
}

/// Test the sized index range
GTEST_TEST(detray_core, sized_index_range) {
  // Index range that contains lower index and number of elements
  using index_t = detail::index_range<unsigned int, detail::sized_index_range>;
  constexpr auto inv_idx{std::numeric_limits<std::uint32_t>::max()};

  static_assert(concepts::interval<index_t>);

  // Check a empty index range (invalid)
  auto ir = index_t{};
  EXPECT_EQ(ir.lower(), inv_idx);
  EXPECT_EQ(ir.upper(), 2u * inv_idx);
  EXPECT_TRUE(ir.is_invalid());

  // Check the values after setting them
  ir.set_lower(5u).set_upper(42u);

  EXPECT_EQ(ir.lower(), 5u);
  EXPECT_EQ(ir.upper(), 42u);
  EXPECT_EQ(detail::get<0>(ir), 5u);
  EXPECT_EQ(detail::get<1>(ir), 42u);
  EXPECT_EQ(ir.size(), 37u);

  // Check invalid link
  EXPECT_FALSE(ir.is_invalid());
  ir.set_lower(inv_idx);
  EXPECT_TRUE(ir.is_invalid());
  ir.set_lower(10u);
  EXPECT_FALSE(ir.is_invalid());

  ir.set_upper(inv_idx);
  EXPECT_TRUE(ir.is_invalid());
  ir.set_upper(20u);
  EXPECT_FALSE(ir.is_invalid());

  EXPECT_EQ(ir.lower(), 10u);
  EXPECT_EQ(ir.upper(), 20u);
  EXPECT_EQ(detail::get<0>(ir), 10u);
  EXPECT_EQ(detail::get<1>(ir), 20u);
  EXPECT_EQ(ir.size(), 10u);

  // Test equality operator
  EXPECT_TRUE(ir == ir);
  index_t other = ir;
  EXPECT_TRUE(ir == other);
  ir.set_lower(13u);
  EXPECT_FALSE(ir == other);

  // Test print operator
  std::clog << ir << std::endl;

  // Test arithmetic operators
  index_t ir2 = ir + index_t{3u, 5u};

  EXPECT_EQ(ir2.lower(), 3u);
  EXPECT_EQ(ir2.upper(), 23);
  EXPECT_EQ(ir2.size(), 20u);
  ir2 = ir2 + 25u;
  EXPECT_EQ(ir2.lower(), 3u);
  EXPECT_EQ(ir2.upper(), 48u);
  EXPECT_EQ(ir2.size(), 45u);

  index_t ir3 = ir - index_t{3u, 15u};

  EXPECT_EQ(ir3.lower(), 13u);
  EXPECT_EQ(ir3.upper(), 18u);
  EXPECT_EQ(ir3.size(), 5u);
  ir3 = ir3 + 2u;
  EXPECT_EQ(ir3.lower(), 13u);
  EXPECT_EQ(ir3.upper(), 20u);
  EXPECT_EQ(ir3.size(), 7u);

  // Comparison operators
  EXPECT_TRUE(ir3 < ir2);
  EXPECT_TRUE(ir3 != ir2);
  EXPECT_TRUE(ir2 > ir3);

  EXPECT_FALSE(ir3 > ir2);
  EXPECT_FALSE(ir3 == ir2);
  EXPECT_FALSE(ir2 < ir3);
}

/// Test the multi index
GTEST_TEST(detray_core, multi_index) {
  // Index range that contains lower index and number of elements
  using index_t = dmulti_index<unsigned int, 3>;
  constexpr auto inv_idx{detail::invalid_value<unsigned int>()};

  // Check a empty index range (invalid)
  auto mi = index_t{};
  EXPECT_EQ(mi[0], inv_idx);
  EXPECT_EQ(mi[1], inv_idx);
  EXPECT_EQ(mi[2], inv_idx);
  EXPECT_EQ(mi.size(), 3u);
  EXPECT_EQ(index_t::size(), 3u);
  EXPECT_TRUE(mi.is_invalid());

  // Check the values after setting them
  mi[0] = 21u;
  mi[1] = 42u;
  mi[2] = 84u;

  EXPECT_EQ(mi[0], 21u);
  EXPECT_EQ(mi[1], 42u);
  EXPECT_EQ(mi[2], 84u);
  EXPECT_EQ(detail::get<0>(mi), 21u);
  EXPECT_EQ(detail::get<1>(mi), 42u);
  EXPECT_EQ(detail::get<2>(mi), 84u);

  // Check invalid link
  EXPECT_FALSE(mi.is_invalid());
  mi[1] = inv_idx;
  EXPECT_TRUE(mi.is_invalid());

  // Test equality operator
  EXPECT_TRUE(mi == mi);
  index_t other = mi;
  EXPECT_TRUE(mi == other);
  mi[1] = 13u;
  EXPECT_FALSE(mi == other);

  // Test print operator
  std::clog << mi << std::endl;
}

/// Test the typed index
GTEST_TEST(detray_core, typed_index) {
  using index_t = dtyped_index<mask_id, unsigned int>;

  // Check a empty index
  auto ti = index_t{};
  EXPECT_EQ(ti.id(), static_cast<mask_id>((1u << 4) - 1u));
  EXPECT_EQ(ti.index(), static_cast<unsigned int>((1u << 28) - 1u));
  EXPECT_TRUE(ti.is_invalid());

  // Check the values after setting them
  ti.set_id(mask_id::e_unmasked).set_index(42u);

  EXPECT_EQ(ti.id(), mask_id::e_unmasked);
  EXPECT_EQ(ti.index(), 42u);

  // Check invalid link
  EXPECT_FALSE(ti.is_invalid());
  ti.set_id(static_cast<index_t::id_type>((1u << 4) - 1u));
  EXPECT_TRUE(ti.is_invalid());
  ti.set_id(mask_id::e_unmasked);
  EXPECT_FALSE(ti.is_invalid());
  ti.set_index((1u << 28) - 1u);
  EXPECT_TRUE(ti.is_invalid());
  ti.set_index(21u);
  EXPECT_FALSE(ti.is_invalid());

  // Test equality operator
  EXPECT_TRUE(ti == ti);
  index_t other = ti;
  EXPECT_TRUE(ti == other);
  ti.set_index(13u);
  EXPECT_FALSE(ti == other);

  // Test print operator
  std::clog << ti << std::endl;

  // Test arithmetic operators
  ti += 5u;

  EXPECT_EQ(ti.id(), mask_id::e_unmasked);
  EXPECT_EQ(ti.index(), 18u);

  ti -= 6u;
  EXPECT_EQ(ti.index(), 12u);

  index_t ti2 = ti + index_t{mask_id::e_unmasked, 6u};
  EXPECT_EQ(ti2.id(), mask_id::e_unmasked);
  EXPECT_EQ(ti2.index(), 18u);
  ++ti2;
  EXPECT_EQ(ti2.index(), 19u);

  ti2 += ti;
  EXPECT_EQ(ti2.id(), mask_id::e_unmasked);
  EXPECT_EQ(ti2.index(), 31u);

  // Comparison operators
  EXPECT_TRUE(ti < ti2);
  EXPECT_TRUE(ti != ti2);
  EXPECT_TRUE(ti2 > ti);

  EXPECT_FALSE(ti > ti2);
  EXPECT_FALSE(ti == ti2);
  EXPECT_FALSE(ti2 < ti);
}

/// Test the typed index using index ranges
GTEST_TEST(detray_core, typed_index_range) {
  using ranged_index_t =
      detail::index_range<unsigned int, detail::sized_index_range,
                          std::uint_least32_t, 0x0fffff00, 0x000000ff>;
  using index_t = dtyped_index<mask_id, ranged_index_t>;

  // Check a empty index
  auto tri = index_t{};
  EXPECT_EQ(tri.id(), static_cast<mask_id>((1u << 4) - 1u));
  EXPECT_EQ(tri.index(), ranged_index_t{});
  EXPECT_TRUE(tri.is_invalid());

  // Check the values after setting them
  constexpr ranged_index_t test_range{10u, 42u};
  tri.set_id(mask_id::e_unmasked).set_index(test_range);

  EXPECT_EQ(tri.id(), mask_id::e_unmasked);
  EXPECT_EQ(tri.index(), test_range);

  // Check invalid link
  EXPECT_FALSE(tri.is_invalid());
  tri.set_id(static_cast<index_t::id_type>((1u << 4) - 1u));
  EXPECT_TRUE(tri.is_invalid());
  tri.set_id(mask_id::e_unmasked);
  EXPECT_FALSE(tri.is_invalid());
  tri.set_index(ranged_index_t{});
  EXPECT_TRUE(tri.is_invalid());
  tri.set_index(ranged_index_t{5u, 42u});
  EXPECT_FALSE(tri.is_invalid());

  // Test equality operator
  EXPECT_TRUE(tri == tri);
  index_t other = tri;
  EXPECT_TRUE(tri == other);
  tri.set_index(ranged_index_t{4u, 13u});
  EXPECT_FALSE(tri == other);

  // Test print operator
  std::clog << tri << std::endl;

  // Test arithmetic operators
  tri += 5u;

  EXPECT_EQ(tri.id(), mask_id::e_unmasked);
  ranged_index_t test_range2{4u, 18u};
  EXPECT_EQ(tri.index(), test_range2);

  tri -= 6u;
  ranged_index_t test_range3{4u, 12u};
  EXPECT_EQ(tri.index(), test_range3);

  index_t tri2 = tri + index_t{mask_id::e_unmasked, ranged_index_t{0u, 6u}};
  ranged_index_t test_range4{0u, 16u};
  EXPECT_EQ(tri2.id(), mask_id::e_unmasked);
  EXPECT_EQ(tri2.index(), test_range4);
  ++tri2;
  test_range4 = test_range4 + 1u;
  EXPECT_EQ(tri2.index(), test_range4);

  tri2 += tri;
  EXPECT_EQ(tri2.id(), mask_id::e_unmasked);
  EXPECT_EQ(tri2.index(), test_range4);

  // Comparison operators
  EXPECT_TRUE(tri < tri2);
  EXPECT_TRUE(tri != tri2);
  EXPECT_TRUE(tri2 > tri);

  EXPECT_FALSE(tri > tri2);
  EXPECT_FALSE(tri == tri2);
  EXPECT_FALSE(tri2 < tri);
}
