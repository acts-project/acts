// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/detail/SortByFloat.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

using namespace Acts;

namespace ActsTests {

namespace {

struct Item {
  std::size_t index{};
  float key{};

  bool operator==(const Item& other) const = default;
};

/// Sorts @p keys with detail::fillAndSortByFloat and checks the result against
/// std::ranges::stable_sort.
void checkAgainstStableSort(const std::vector<float>& keys) {
  const auto itemAt = [&](std::size_t k) { return Item{k, keys[k]}; };
  const auto keyOf = [](const Item& item) { return item.key; };

  std::vector<Item> expected;
  for (std::size_t k = 0; k < keys.size(); ++k) {
    expected.push_back(itemAt(k));
  }
  std::ranges::stable_sort(expected, {}, keyOf);

  std::vector<Item> sorted;
  std::vector<Item> scratch;
  detail::fillAndSortByFloat(sorted, scratch, keys.size(), itemAt, keyOf);

  BOOST_CHECK(sorted == expected);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(SortByFloat_MatchesStableSort) {
  std::mt19937 rng(42);
  std::uniform_real_distribution<float> uniform(-5.f, 5.f);
  // Sizes on both sides of the rank sort, radix sort and fallback limits
  for (std::size_t n : {0, 1, 2, 3, 39, 40, 41, 100, 1023, 1024, 1025, 3000}) {
    std::vector<float> keys(n);
    for (float& key : keys) {
      key = uniform(rng);
    }
    checkAgainstStableSort(keys);
  }
}

BOOST_AUTO_TEST_CASE(SortByFloat_KeepsOrderOfEqualKeys) {
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> fewValues(-3, 3);
  for (std::size_t n : {10, 40, 41, 500, 2000}) {
    std::vector<float> keys(n);
    for (float& key : keys) {
      key = static_cast<float>(fewValues(rng));
    }
    checkAgainstStableSort(keys);
  }
  checkAgainstStableSort(std::vector<float>(30, 1.f));
  checkAgainstStableSort(std::vector<float>(300, 1.f));
}

BOOST_AUTO_TEST_CASE(SortByFloat_UnusualKeys) {
  const float inf = std::numeric_limits<float>::infinity();
  const float largest = std::numeric_limits<float>::max();
  const float smallest = std::numeric_limits<float>::denorm_min();
  for (std::size_t n : {20, 200}) {
    std::vector<float> keys(n);
    for (std::size_t k = 0; k < n; ++k) {
      keys[k] = static_cast<float>(n - k);
    }
    // A key that is not finite
    std::vector<float> withInf = keys;
    withInf[n / 2] = inf;
    withInf[n / 3] = -inf;
    checkAgainstStableSort(withInf);
    // Keys too large to scale into buckets
    std::vector<float> large = keys;
    large[n / 2] = largest;
    large[n / 3] = -largest;
    checkAgainstStableSort(large);
    // Keys too close together to scale into buckets
    std::vector<float> close(n);
    for (std::size_t k = 0; k < n; ++k) {
      close[k] = (k % 2 == 0) ? 0.f : smallest;
    }
    checkAgainstStableSort(close);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
