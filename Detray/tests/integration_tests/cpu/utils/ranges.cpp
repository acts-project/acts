// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// detray core
#include "detray/utils/ranges.hpp"

#include "detray/definitions/containers.hpp"

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// Integration test for enumeration of a subrange
GTEST_TEST(detray_utils, ranges_subrange_iota) {
  std::array<std::size_t, 2> seq{1u, 10u};
  std::array<std::size_t, 2> interval{3u, 7u};

  auto iota_sr = detray::ranges::subrange(detray::views::iota(seq), interval);

  // Check iterator category
  static_assert(detray::ranges::forward_range<decltype(iota_sr)>);
  static_assert(detray::ranges::bidirectional_range<decltype(iota_sr)>);
  static_assert(detray::ranges::random_access_range<decltype(iota_sr)>);

  std::size_t i{4u};
  for (const auto v : iota_sr) {
    ASSERT_EQ(i++, v);
  }

  // Check range composition with pipe operator
  auto iota_sr_comp = detray::views::iota(seq) |
                      detray::ranges::subrange(interval[0], interval[1]);
  i = 4u;
  for (const auto v : iota_sr_comp) {
    ASSERT_EQ(i++, v);
  }
}

// Integration test for enumeration of a subrange
GTEST_TEST(detray_utils, ranges_enumerated_subrange) {
  struct uint_holder {
    unsigned int ui{0u};
  };

  dvector<uint_holder> seq = {{0u}, {1u}, {2u}, {3u}, {4u}, {5u}};

  std::size_t begin{1u};
  std::size_t end{4u};
  std::array<std::size_t, 2> interval{begin, end};

  auto enum_sr =
      detray::views::enumerate(detray::ranges::subrange(seq, interval));

  // Check iterator category
  static_assert(detray::ranges::random_access_range<decltype(enum_sr)>);

  dvector<unsigned int> expected{};
  for (const auto [i, v] : enum_sr) {
    ASSERT_EQ(i, v.ui - 1u);
    expected.push_back(v.ui);
  }

  // Check range composition with pipe operator
  auto enum_sr_comp = seq | detray::ranges::subrange(interval[0], interval[1]) |
                      detray::views::enumerate(2u);

  for (const auto [i, v] : enum_sr_comp) {
    // Enumeration starts at 'two'
    ASSERT_EQ(expected[i - 2u], v.ui);
  }
}

// Integration test for the picking of indexed elements from another range
GTEST_TEST(detray_utils, ranges_pick_static_joined_sequence) {
  darray<dindex, 2> interval_1 = {2u, 4u};
  darray<dindex, 2> interval_2 = {7u, 9u};

  // The indices of the iota elements to be picked
  std::vector<dindex> reference = {2u, 3u, 7u, 8u};
  std::vector<dindex> check = {};

  struct uint_holder {
    unsigned int ui{0u};
  };

  dvector<uint_holder> seq = {{0u}, {1u}, {2u}, {3u}, {4u},
                              {5u}, {6u}, {7u}, {8u}};

  auto indices = detray::views::static_join(detray::views::iota(interval_1),
                                            detray::views::iota(interval_2));
  auto selected = detray::views::pick(seq, indices);

  // Check iterator category
  static_assert(detray::ranges::forward_range<decltype(indices)>);
  static_assert(detray::ranges::bidirectional_range<decltype(indices)>);
  static_assert(detray::ranges::random_access_range<decltype(indices)>);
  static_assert(detray::ranges::forward_range<decltype(selected)>);
  static_assert(detray::ranges::bidirectional_range<decltype(selected)>);
  static_assert(detray::ranges::random_access_range<decltype(selected)>);

  // Test inherited member functions
  const auto [i, v] = selected[2];
  ASSERT_EQ(i, 7u);
  ASSERT_EQ(v.ui, 7u);
  ASSERT_EQ(selected.size(), 4u);
  const auto [i_front, v_front] = selected.front();
  ASSERT_EQ(i_front, 2u);
  ASSERT_EQ(v_front.ui, 2u);
  const auto [i_back, v_back] = selected.back();
  ASSERT_EQ(i_back, 8u);
  ASSERT_EQ(v_back.ui, 8u);

  for (auto [j, w] : selected) {
    ASSERT_TRUE(j == w.ui);
    check.push_back(w.ui);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);

  // Check range composition with pipe operator
  auto selected_comp = seq | detray::views::pick(indices);

  check.clear();
  for (auto [j, w] : selected_comp) {
    ASSERT_TRUE(j == w.ui);
    check.push_back(w.ui);
  }
  ASSERT_EQ(check.size(), reference.size());
  ASSERT_EQ(check, reference);
}
