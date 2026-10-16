// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/detail/SortByFloat.hpp"

#include <algorithm>
#include <cmath>
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
};

const auto keyOf = [](const Item& item) { return item.key; };

std::vector<Item> itemsOf(const std::vector<float>& keys) {
  std::vector<Item> items;
  for (std::size_t k = 0; k < keys.size(); ++k) {
    items.push_back(Item{k, keys[k]});
  }
  return items;
}

/// Checks that @p sorted holds the items of @p keys in the order
/// std::ranges::stable_sort gives, keys bit for bit (the sign of zero
/// included).
void checkSorted(const std::vector<float>& keys,
                 const std::vector<Item>& sorted) {
  std::vector<Item> expected = itemsOf(keys);
  std::ranges::stable_sort(expected, {}, keyOf);

  std::vector<std::size_t> expectedIndices;
  std::vector<std::size_t> sortedIndices;
  for (const Item& item : expected) {
    expectedIndices.push_back(item.index);
  }
  for (const Item& item : sorted) {
    sortedIndices.push_back(item.index);
    BOOST_CHECK_EQUAL(std::signbit(item.key), std::signbit(keys[item.index]));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(sortedIndices.begin(), sortedIndices.end(),
                                expectedIndices.begin(), expectedIndices.end());
}

/// Sorts @p keys with detail::fillAndSortByFloat, reusing @p scratch, and
/// checks the result.
void checkFillAndSort(const std::vector<float>& keys,
                      std::vector<Item>& scratch) {
  const auto itemAt = [&](std::size_t k) { return Item{k, keys[k]}; };
  std::vector<Item> sorted;
  detail::fillAndSortByFloat(sorted, scratch, keys.size(), itemAt, keyOf);
  checkSorted(keys, sorted);
}

std::vector<float> uniformKeys(std::size_t n, std::mt19937& rng) {
  std::uniform_real_distribution<float> uniform(-5.f, 5.f);
  std::vector<float> keys(n);
  for (float& key : keys) {
    key = uniform(rng);
  }
  return keys;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(SortByFloat_MatchesStableSort) {
  std::mt19937 rng(42);
  // One scratch for every call, of whatever size the last call left it
  std::vector<Item> scratch;
  // Sizes on both sides of the insertion, rank and bucket sort limits
  for (std::size_t n : {0,  1,  2,  3,  4,   5,    6,    7,    8,    9,
                        10, 39, 40, 41, 100, 1023, 1024, 1025, 3000, 41}) {
    checkFillAndSort(uniformKeys(n, rng), scratch);
  }
}

BOOST_AUTO_TEST_CASE(SortByFloat_EachSortMatchesStableSort) {
  std::mt19937 rng(3);
  std::vector<Item> scratch;
  for (std::size_t n : {0, 1, 2, 7, 8, 9, 39, 40}) {
    const std::vector<float> keys = uniformKeys(n, rng);
    std::vector<Item> items = itemsOf(keys);
    detail::insertionSortByFloat(items, keyOf);
    checkSorted(keys, items);
    items = itemsOf(keys);
    detail::rankSortByFloat<40>(items, keyOf);
    checkSorted(keys, items);
  }
  for (std::size_t n : {2, 41, 100, 1024, 1025}) {
    const std::vector<float> keys = uniformKeys(n, rng);
    std::vector<Item> items = itemsOf(keys);
    const auto [minKey, maxKey] = std::ranges::minmax(keys);
    detail::bucketSortByFloat(items, scratch, keyOf, minKey, maxKey);
    checkSorted(keys, items);
  }
}

BOOST_AUTO_TEST_CASE(SortByFloat_KeepsOrderOfEqualKeys) {
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> fewValues(-3, 3);
  std::vector<Item> scratch;
  for (std::size_t n : {5, 10, 40, 41, 500, 2000}) {
    std::vector<float> keys(n);
    for (float& key : keys) {
      key = static_cast<float>(fewValues(rng));
    }
    checkFillAndSort(keys, scratch);
  }
  checkFillAndSort(std::vector<float>(6, 1.f), scratch);
  checkFillAndSort(std::vector<float>(30, 1.f), scratch);
  checkFillAndSort(std::vector<float>(300, 1.f), scratch);
}

BOOST_AUTO_TEST_CASE(SortByFloat_SignedZeros) {
  std::vector<Item> scratch;
  for (std::size_t n : {6, 20, 200}) {
    std::vector<float> keys(n);
    for (std::size_t k = 0; k < n; ++k) {
      keys[k] = (k % 3 == 0) ? -0.f : (k % 3 == 1 ? 0.f : 1.f - 0.01f * k);
    }
    checkFillAndSort(keys, scratch);
  }
}

BOOST_AUTO_TEST_CASE(SortByFloat_UnusualKeys) {
  const float inf = std::numeric_limits<float>::infinity();
  const float largest = std::numeric_limits<float>::max();
  const float smallest = std::numeric_limits<float>::denorm_min();
  std::mt19937 rng(11);
  std::vector<Item> scratch;
  for (std::size_t n : {6, 20, 200, 1000}) {
    std::vector<float> keys(n);
    for (std::size_t k = 0; k < n; ++k) {
      keys[k] = static_cast<float>(n - k);
    }
    // A key that is not finite
    std::vector<float> withInf = keys;
    withInf[n / 2] = inf;
    withInf[n / 3] = -inf;
    checkFillAndSort(withInf, scratch);
    // Keys too large to scale into buckets
    std::vector<float> large = keys;
    large[n / 2] = largest;
    large[n / 3] = -largest;
    checkFillAndSort(large, scratch);
    // Keys too close together to scale into buckets
    std::vector<float> close(n);
    for (std::size_t k = 0; k < n; ++k) {
      close[k] = (k % 2 == 0) ? 0.f : smallest;
    }
    checkFillAndSort(close, scratch);
    // One outlier: almost every item in one bucket
    std::uniform_real_distribution<float> narrow(-0.005f, 0.005f);
    std::vector<float> outlier(n);
    for (float& key : outlier) {
      key = narrow(rng);
    }
    outlier[n / 2] = 27.f;
    checkFillAndSort(outlier, scratch);
  }
}

BOOST_AUTO_TEST_CASE(SortByFloat_SortByCotThetaOverloads) {
  std::mt19937 rng(5);
  std::uniform_real_distribution<float> cotTheta(-4.f, 4.f);
  for (std::size_t n : {0, 5, 30, 300}) {
    DoubletsForMiddleSp doublets;
    for (std::size_t k = 0; k < n; ++k) {
      doublets.emplace_back(static_cast<SpacePointIndex>(k), cotTheta(rng), 1.f,
                            1.f, 0.f, 0.f, 0.f, 0.f);
    }
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> allocating;
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> withScratch;
    std::vector<DoubletsForMiddleSp::IndexAndCotTheta> scratch;
    doublets.sortByCotTheta({0, doublets.size()}, allocating);
    doublets.sortByCotTheta({0, doublets.size()}, withScratch, scratch);
    BOOST_REQUIRE_EQUAL(allocating.size(), n);
    BOOST_REQUIRE_EQUAL(withScratch.size(), n);
    for (std::size_t k = 0; k < n; ++k) {
      BOOST_CHECK_EQUAL(allocating[k].index, withScratch[k].index);
      if (k > 0) {
        BOOST_CHECK(!(withScratch[k].cotTheta < withScratch[k - 1].cotTheta));
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
