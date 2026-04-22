// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray core
#include "detray/utils/grid/populators.hpp"

#include "detray/definitions/indexing.hpp"
#include "detray/utils/grid/grid_bins.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <climits>

using namespace detray;

namespace {

constexpr dindex inf{std::numeric_limits<dindex>::max()};

// Create some bin data for non-owning grid
template <typename bin_t>
struct bin_content_sequence {
  using entry_t = typename bin_t::entry_type;
  entry_t entry{0};

  auto operator()() {
    entry += entry_t{1};
    return bin_t{}.init(entry);
  }
};

/// Test single entry bin content element by element
template <typename bin_content_t, typename content_t>
void test_content(const bin_content_t& bin_content, const content_t& expected) {
  unsigned int i{0u};
  for (const auto& elem : bin_content) {
    ASSERT_EQ(elem, expected[i++]);
  }
}

}  // anonymous namespace

/// Replace populator
GTEST_TEST(detray_grid, replace_populator) {
  detray::replace replacer{};

  // Create some bin data
  dvector<bins::single<dindex>> bin_data{};
  bin_data.resize(50);
  std::ranges::generate_n(bin_data.begin(), 50,
                          bin_content_sequence<bins::single<dindex>>());

  // Check test setup
  EXPECT_EQ(bin_data[0].value(), 1u);
  EXPECT_EQ(bin_data[1].value(), 2u);
  EXPECT_EQ(bin_data[2].value(), 3u);
  EXPECT_EQ(bin_data[3].value(), 4u);
  EXPECT_EQ(bin_data[19].value(), 20u);
  EXPECT_EQ(bin_data[20].value(), 21u);
  EXPECT_EQ(bin_data[21].value(), 22u);
  EXPECT_EQ(bin_data[41].value(), 42u);
  EXPECT_EQ(bin_data[42].value(), 43u);
  EXPECT_EQ(bin_data[43].value(), 44u);
  EXPECT_EQ(bin_data[48].value(), 49u);
  EXPECT_EQ(bin_data[49].value(), 50u);

  // Replace some bin entries
  replacer(bin_data[0], 500u);
  replacer(bin_data[2], 5u);
  replacer(bin_data[20], 50u);
  replacer(bin_data[42], 55u);
  replacer(bin_data[49], 6u);

  EXPECT_EQ(bin_data[0].value(), 500u);
  EXPECT_EQ(bin_data[1].value(), 2u);
  EXPECT_EQ(bin_data[2].value(), 5u);
  EXPECT_EQ(bin_data[3].value(), 4u);
  EXPECT_EQ(bin_data[19].value(), 20u);
  EXPECT_EQ(bin_data[20].value(), 50u);
  EXPECT_EQ(bin_data[21].value(), 22u);
  EXPECT_EQ(bin_data[41].value(), 42u);
  EXPECT_EQ(bin_data[42].value(), 55u);
  EXPECT_EQ(bin_data[43].value(), 44u);
  EXPECT_EQ(bin_data[48].value(), 49u);
  EXPECT_EQ(bin_data[49].value(), 6u);
}

/// Complete populator
GTEST_TEST(detray_grid, complete_populator) {
  detray::complete<> completer{};

  // Create some bin data
  dvector<bins::static_array<dindex, 4>> bin_data{};
  bin_data.resize(50);
  std::ranges::generate_n(
      bin_data.begin(), 50,
      bin_content_sequence<bins::static_array<dindex, 4>>());

  // Check test setup
  std::array<dindex, 4> stored = {1u, inf, inf, inf};
  test_content(bin_data[0], stored);
  stored = {2u, inf, inf, inf};
  test_content(bin_data[1], stored);
  stored = {3u, inf, inf, inf};
  test_content(bin_data[2], stored);
  stored = {4u, inf, inf, inf};
  test_content(bin_data[3], stored);
  stored = {20u, inf, inf, inf};
  test_content(bin_data[19], stored);
  stored = {21u, inf, inf, inf};
  test_content(bin_data[20], stored);
  stored = {22u, inf, inf, inf};
  test_content(bin_data[21], stored);
  stored = {42u, inf, inf, inf};
  test_content(bin_data[41], stored);
  stored = {43u, inf, inf, inf};
  test_content(bin_data[42], stored);
  stored = {44u, inf, inf, inf};
  test_content(bin_data[43], stored);
  stored = {49u, inf, inf, inf};
  test_content(bin_data[48], stored);
  stored = {50u, inf, inf, inf};
  test_content(bin_data[49], stored);

  // Fill the bins with some bin entry
  completer(bin_data[0], 500u);
  completer(bin_data[2], 5u);
  completer(bin_data[20], 50u);
  completer(bin_data[42], 55u);
  completer(bin_data[49], 6u);

  stored = {1u, 500u, 500u, 500u};
  test_content(bin_data[0], stored);
  stored = {2u, inf, inf, inf};
  test_content(bin_data[1], stored);
  stored = {3u, 5u, 5u, 5u};
  test_content(bin_data[2], stored);
  stored = {4u, inf, inf, inf};
  test_content(bin_data[3], stored);
  stored = {20u, inf, inf, inf};
  test_content(bin_data[19], stored);
  stored = {21u, 50u, 50u, 50u};
  test_content(bin_data[20], stored);
  stored = {22u, inf, inf, inf};
  test_content(bin_data[21], stored);
  stored = {42u, inf, inf, inf};
  test_content(bin_data[41], stored);
  stored = {43u, 55u, 55u, 55u};
  test_content(bin_data[42], stored);
  stored = {44u, inf, inf, inf};
  test_content(bin_data[43], stored);
  stored = {49u, inf, inf, inf};
  test_content(bin_data[48], stored);
  stored = {50u, 6u, 6u, 6u};
  test_content(bin_data[49], stored);
}

/// Regular attach populator
GTEST_TEST(detray_grid, regular_attach_populator) {
  detray::attach<> attacher{};

  // Create some bin data
  dvector<bins::static_array<dindex, 4>> bin_data{};
  bin_data.resize(50);
  std::ranges::generate_n(
      bin_data.begin(), 50,
      bin_content_sequence<bins::static_array<dindex, 4>>());

  // Check test setup
  std::array<dindex, 4> stored = {1u, inf, inf, inf};
  test_content(bin_data[0], stored);
  stored = {2u, inf, inf, inf};
  test_content(bin_data[1], stored);
  stored = {3u, inf, inf, inf};
  test_content(bin_data[2], stored);
  stored = {4u, inf, inf, inf};
  test_content(bin_data[3], stored);
  stored = {20u, inf, inf, inf};
  test_content(bin_data[19], stored);
  stored = {21u, inf, inf, inf};
  test_content(bin_data[20], stored);
  stored = {22u, inf, inf, inf};
  test_content(bin_data[21], stored);
  stored = {42u, inf, inf, inf};
  test_content(bin_data[41], stored);
  stored = {43u, inf, inf, inf};
  test_content(bin_data[42], stored);
  stored = {44u, inf, inf, inf};
  test_content(bin_data[43], stored);
  stored = {49u, inf, inf, inf};
  test_content(bin_data[48], stored);
  stored = {50u, inf, inf, inf};
  test_content(bin_data[49], stored);

  // Fill the bins with some bin entry
  attacher(bin_data[0], 500u);
  attacher(bin_data[2], 5u);
  attacher(bin_data[2], 79u);
  attacher(bin_data[20], 50u);
  attacher(bin_data[42], 55u);
  attacher(bin_data[49], 6u);
  attacher(bin_data[49], 7u);
  attacher(bin_data[49], 8u);

  stored = {1u, 500u, inf, inf};
  test_content(bin_data[0], stored);
  stored = {2u, inf, inf, inf};
  test_content(bin_data[1], stored);
  stored = {3u, 5u, 79u, inf};
  test_content(bin_data[2], stored);
  stored = {4u, inf, inf, inf};
  test_content(bin_data[3], stored);
  stored = {20u, inf, inf, inf};
  test_content(bin_data[19], stored);
  stored = {21u, 50u, inf, inf};
  test_content(bin_data[20], stored);
  stored = {22u, inf, inf, inf};
  test_content(bin_data[21], stored);
  stored = {42u, inf, inf, inf};
  test_content(bin_data[41], stored);
  stored = {43u, 55u, inf, inf};
  test_content(bin_data[42], stored);
  stored = {44u, inf, inf, inf};
  test_content(bin_data[43], stored);
  stored = {49u, inf, inf, inf};
  test_content(bin_data[48], stored);
  stored = {50u, 6u, 7u, 8u};
  test_content(bin_data[49], stored);
}
