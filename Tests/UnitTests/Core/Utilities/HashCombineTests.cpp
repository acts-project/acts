// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/HashCombine.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_set>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(hashMix_nonIdentity) {
  // hashMix should not be the identity for small integers.
  // Note: fmix64(0) == 0 is a known fixed point of the murmur3 finalizer.
  for (std::size_t i = 1; i < 100; ++i) {
    BOOST_CHECK_NE(hashMix(i), i);
  }
}

BOOST_AUTO_TEST_CASE(hashMix_distinct) {
  // consecutive inputs should produce distinct outputs
  std::unordered_set<std::size_t> seen;
  for (std::size_t i = 0; i < 10000; ++i) {
    auto [it, inserted] = seen.insert(hashMix(i));
    BOOST_CHECK(inserted);
  }
}

BOOST_AUTO_TEST_CASE(hashMixAndCombine_basic) {
  // same inputs produce same hash
  std::size_t a = hashMixAndCombine(1, 2, 3);
  std::size_t b = hashMixAndCombine(1, 2, 3);
  BOOST_CHECK_EQUAL(a, b);

  // different inputs produce different hashes
  std::size_t c = hashMixAndCombine(1, 2, 4);
  BOOST_CHECK_NE(a, c);
}

BOOST_AUTO_TEST_CASE(hashMixAndCombine_orderMatters) {
  std::size_t a = hashMixAndCombine(1, 2);
  std::size_t b = hashMixAndCombine(2, 1);
  BOOST_CHECK_NE(a, b);
}

BOOST_AUTO_TEST_CASE(hashMixAndCombine_consecutiveCollisionRate) {
  // Simulate barcode-like inputs: 5 small integer fields, many with zeros.
  // Hash into a 4096-bucket table and check collision quality.
  constexpr std::size_t nBuckets = 4096;

  std::unordered_set<std::size_t> hashes;
  std::array<std::size_t, nBuckets> buckets{};

  // Generate barcode-like inputs: particle ID 0..99999, across several
  // vertex primaries, with fixed other fields.
  std::size_t nInputs = 0;
  for (std::uint32_t vp = 0; vp < 10; ++vp) {
    for (std::uint32_t p = 0; p < 100000; ++p) {
      std::size_t h = hashMixAndCombine(vp, std::uint32_t{0}, p,
                                        std::uint32_t{0}, std::uint32_t{0});
      hashes.insert(h);
      buckets[h % nBuckets]++;
      ++nInputs;
    }
  }

  // No full collisions for 1M distinct inputs
  BOOST_CHECK_EQUAL(hashes.size(), nInputs);

  // Check that the distribution across buckets is reasonable.
  // With 1M items in 4096 buckets (~244 per bucket on average), a good hash
  // should keep the max bucket well below 2x the average.
  std::size_t maxBucket = *std::max_element(buckets.begin(), buckets.end());
  double avgBucket = static_cast<double>(nInputs) / nBuckets;
  double ratio = static_cast<double>(maxBucket) / avgBucket;
  BOOST_CHECK_LT(ratio, 2.0);
}

BOOST_AUTO_TEST_CASE(hashMixAndCombine_multiFieldCollisionRate) {
  // All fields vary in small ranges (realistic barcode space)
  std::unordered_set<std::size_t> hashes;
  std::size_t nInputs = 0;

  // 20 * 10 * 500 * 5 * 5 = 2'500'000 inputs
  for (std::uint32_t vp = 0; vp < 20; ++vp) {
    for (std::uint32_t vs = 0; vs < 10; ++vs) {
      for (std::uint32_t p = 0; p < 500; ++p) {
        for (std::uint32_t g = 0; g < 5; ++g) {
          for (std::uint32_t sp = 0; sp < 5; ++sp) {
            std::size_t h = hashMixAndCombine(vp, vs, p, g, sp);
            hashes.insert(h);
            ++nInputs;
          }
        }
      }
    }
  }

  // All 2.5M inputs should produce unique hashes
  BOOST_CHECK_EQUAL(hashes.size(), nInputs);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
