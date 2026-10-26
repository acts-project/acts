// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/DoubletSeedFinder.hpp"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

using namespace Acts;

namespace ActsTests {

namespace {

using IndexAndCotTheta = DoubletsForMiddleSp::IndexAndCotTheta;

DoubletsForMiddleSp doubletsOf(const std::vector<float>& cotThetas) {
  DoubletsForMiddleSp doublets;
  for (std::size_t k = 0; k < cotThetas.size(); ++k) {
    doublets.emplace_back(static_cast<SpacePointIndex>(k), cotThetas[k], 1.f,
                          1.f, 0.f, 0.f, 0.f, 0.f);
  }
  return doublets;
}

/// Checks that @p sorted holds the doublets of @p range in the order
/// std::ranges::stable_sort gives, with the sign of each cotTheta kept (which
/// tells -0 from +0).
void checkSorted(const std::vector<float>& cotThetas,
                 const DoubletsForMiddleSp::IndexRange& range,
                 const std::vector<IndexAndCotTheta>& sorted) {
  std::vector<IndexAndCotTheta> expected;
  for (auto i = range.first; i < range.second; ++i) {
    expected.push_back({i, cotThetas[i]});
  }
  std::ranges::stable_sort(expected, {}, &IndexAndCotTheta::cotTheta);

  std::vector<DoubletsForMiddleSp::Index> expectedIndices;
  std::vector<DoubletsForMiddleSp::Index> sortedIndices;
  std::vector<bool> expectedSigns;
  std::vector<bool> sortedSigns;
  for (const IndexAndCotTheta& item : expected) {
    expectedIndices.push_back(item.index);
    expectedSigns.push_back(std::signbit(item.cotTheta));
  }
  for (const IndexAndCotTheta& item : sorted) {
    sortedIndices.push_back(item.index);
    sortedSigns.push_back(std::signbit(item.cotTheta));
  }

  BOOST_CHECK_EQUAL_COLLECTIONS(sortedIndices.begin(), sortedIndices.end(),
                                expectedIndices.begin(), expectedIndices.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(sortedSigns.begin(), sortedSigns.end(),
                                expectedSigns.begin(), expectedSigns.end());
}

/// Sorts @p range of @p doublets into @p sorted, and checks that no floating
/// point exception was raised: it would trap where they are enabled.
void sortWithoutFpe(const DoubletsForMiddleSp& doublets,
                    const DoubletsForMiddleSp::IndexRange& range,
                    std::vector<IndexAndCotTheta>& sorted,
                    DoubletsForMiddleSp::SortStorage* storage) {
  std::feclearexcept(FE_ALL_EXCEPT);
  if (storage != nullptr) {
    doublets.sortByCotTheta(range, sorted, *storage);
  } else {
    doublets.sortByCotTheta(range, sorted);
  }
  BOOST_CHECK(!std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
}

/// Sorts all of @p cotThetas, reusing @p storage, and checks the result.
void checkSort(const std::vector<float>& cotThetas,
               DoubletsForMiddleSp::SortStorage& storage) {
  const DoubletsForMiddleSp doublets = doubletsOf(cotThetas);
  const DoubletsForMiddleSp::IndexRange range{0, doublets.size()};
  std::vector<IndexAndCotTheta> sorted;
  sortWithoutFpe(doublets, range, sorted, &storage);
  checkSorted(cotThetas, range, sorted);
}

std::vector<float> uniformCotThetas(std::size_t n, std::mt19937& rng) {
  std::uniform_real_distribution<float> uniform(-5.f, 5.f);
  std::vector<float> cotThetas(n);
  for (float& cotTheta : cotThetas) {
    cotTheta = uniform(rng);
  }
  return cotThetas;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(SortByCotThetaMatchesStableSort) {
  std::mt19937 rng(42);
  // One storage for every call, of whatever size the last call left it
  DoubletsForMiddleSp::SortStorage storage;
  // Sizes on both sides of where the sort changes
  for (std::size_t n : {0,  1,  2,  3,  4,  5,  6,  7,  8,   9,    10,   15,
                        16, 17, 31, 32, 33, 40, 64, 65, 100, 1024, 3000, 33}) {
    checkSort(uniformCotThetas(n, rng), storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaKeepsOrderOfEqualCotTheta) {
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> fewValues(-3, 3);
  DoubletsForMiddleSp::SortStorage storage;
  for (std::size_t n : {5, 10, 40, 41, 500, 2000}) {
    std::vector<float> cotThetas(n);
    for (float& cotTheta : cotThetas) {
      cotTheta = static_cast<float>(fewValues(rng));
    }
    checkSort(cotThetas, storage);
  }

  // All equal
  for (std::size_t n : {6, 30, 300}) {
    checkSort(std::vector<float>(n, 1.f), storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaKeepsSignedZeros) {
  DoubletsForMiddleSp::SortStorage storage;
  for (std::size_t n : {6, 20, 200}) {
    std::vector<float> cotThetas(n);
    for (std::size_t k = 0; k < n; ++k) {
      if (k % 3 == 0) {
        cotThetas[k] = -0.f;
      } else if (k % 3 == 1) {
        cotThetas[k] = 0.f;
      } else {
        cotThetas[k] = 1.f - 0.01f * static_cast<float>(k);
      }
    }
    checkSort(cotThetas, storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaFallsBackForInfiniteCotTheta) {
  const float inf = std::numeric_limits<float>::infinity();
  DoubletsForMiddleSp::SortStorage storage;
  for (std::size_t n : {6, 20, 200, 1000}) {
    std::vector<float> cotThetas(n);
    for (std::size_t k = 0; k < n; ++k) {
      cotThetas[k] = static_cast<float>(n - k);
    }
    cotThetas[n / 2] = inf;
    cotThetas[n / 3] = -inf;
    checkSort(cotThetas, storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaFallsBackForTinyRange) {
  // Too close together to spread over buckets
  const float smallest = std::numeric_limits<float>::denorm_min();
  DoubletsForMiddleSp::SortStorage storage;
  for (std::size_t n : {6, 20, 200, 1000}) {
    std::vector<float> cotThetas(n);
    for (std::size_t k = 0; k < n; ++k) {
      cotThetas[k] = (k % 2 == 0) ? 0.f : smallest;
    }
    checkSort(cotThetas, storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaHandlesOneOutlier) {
  // Almost every doublet in one bucket
  std::mt19937 rng(11);
  std::uniform_real_distribution<float> narrow(-0.005f, 0.005f);
  DoubletsForMiddleSp::SortStorage storage;
  for (std::size_t n : {6, 20, 200, 1000}) {
    std::vector<float> cotThetas(n);
    for (float& cotTheta : cotThetas) {
      cotTheta = narrow(rng);
    }
    cotThetas[n / 2] = 27.f;
    checkSort(cotThetas, storage);
  }
}

BOOST_AUTO_TEST_CASE(SortByCotThetaSortsRangesWithBothOverloads) {
  std::mt19937 rng(5);
  for (std::size_t n : {0, 5, 30, 300}) {
    const std::vector<float> cotThetas = uniformCotThetas(n + 20, rng);
    const DoubletsForMiddleSp doublets = doubletsOf(cotThetas);

    // A range that does not start at the first doublet
    const DoubletsForMiddleSp::IndexRange range{
        10, static_cast<DoubletsForMiddleSp::Index>(10 + n)};

    std::vector<IndexAndCotTheta> allocating;
    sortWithoutFpe(doublets, range, allocating, nullptr);
    checkSorted(cotThetas, range, allocating);

    std::vector<IndexAndCotTheta> withStorage;
    DoubletsForMiddleSp::SortStorage storage;
    sortWithoutFpe(doublets, range, withStorage, &storage);
    checkSorted(cotThetas, range, withStorage);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
