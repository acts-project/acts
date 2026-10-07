// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/TrackFinding/SeedDeduplicator.hpp"

#include <limits>
#include <numeric>
#include <span>
#include <vector>

using namespace Acts;

namespace ActsTests {

using Key = SeedDeduplicator::Key;

/// Find the smallest number of seed keys that one track must hold to make a
/// seed of `nKeys` keys a duplicate
Key requiredSharedKeys(const SeedDeduplicator::Config& config, Key nKeys) {
  std::vector<Key> seed(nKeys);
  std::iota(seed.begin(), seed.end(), Key{0});
  for (Key nShared = 1; nShared <= nKeys; ++nShared) {
    SeedDeduplicator dedup(config);
    dedup.reset(nKeys);
    dedup.addTrack(std::span(seed).first(nShared));
    if (dedup.isDuplicate(seed)) {
      return nShared;
    }
  }
  return nKeys + 1;
}

BOOST_AUTO_TEST_SUITE(TrackFindingSuite)

BOOST_AUTO_TEST_CASE(SeedDeduplicatorRequiredSharedScore) {
  // all keys, the default
  BOOST_CHECK_EQUAL(requiredSharedKeys({}, 1), 1u);
  BOOST_CHECK_EQUAL(requiredSharedKeys({}, 3), 3u);

  // N - 1, never below 1
  const SeedDeduplicator::Config relative{.minScore = 0, .maxMissingScore = 1};
  BOOST_CHECK_EQUAL(requiredSharedKeys(relative, 1), 1u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(relative, 3), 2u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(relative, 6), 5u);

  // max(3, N - 2), never above N
  const SeedDeduplicator::Config both{.minScore = 3, .maxMissingScore = 2};
  BOOST_CHECK_EQUAL(requiredSharedKeys(both, 2), 2u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(both, 3), 3u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(both, 5), 3u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(both, 6), 4u);
  BOOST_CHECK_EQUAL(requiredSharedKeys(both, 8), 6u);
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorAllKeys) {
  SeedDeduplicator dedup({});
  dedup.reset(10);

  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 2}));

  dedup.addTrack(std::vector<Key>{0, 1, 2, 3});
  BOOST_CHECK_EQUAL(dedup.nTracks(), 1u);

  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1, 2}));
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{3, 1, 0}));
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 4}));
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{}));
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorTracksDoNotAdd) {
  SeedDeduplicator dedup({.minScore = 0, .maxMissingScore = 1});
  dedup.reset(10);

  // each track holds one key of the seed, the union holds all three
  dedup.addTrack(std::vector<Key>{0, 5});
  dedup.addTrack(std::vector<Key>{1, 6});
  dedup.addTrack(std::vector<Key>{2, 7});
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 2}));

  // a track with two of the three keys
  dedup.addTrack(std::vector<Key>{1, 2, 8});
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1, 2}));
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorSharedKeys) {
  SeedDeduplicator dedup({});
  dedup.reset(10);

  // both tracks hold key 1
  dedup.addTrack(std::vector<Key>{0, 1, 2});
  dedup.addTrack(std::vector<Key>{1, 3, 4});
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1, 2}));
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{1, 3, 4}));
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 3}));
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorRepeatedTrackKey) {
  SeedDeduplicator dedup({});
  dedup.reset(10);

  // a key that occurs twice on a track counts once
  dedup.addTrack(std::vector<Key>{0, 0, 1});
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 2}));
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1}));
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorWeights) {
  // pixel weight 2, strip weight 1, at least 5 shared constraints
  SeedDeduplicator dedup(
      {.minScore = 5,
       .maxMissingScore = std::numeric_limits<SeedDeduplicator::Score>::max()});
  dedup.reset(10);
  for (Key key : {0u, 1u, 2u}) {
    dedup.setWeight(key, 2);
  }
  BOOST_CHECK_EQUAL(dedup.weight(0), 2u);
  BOOST_CHECK_EQUAL(dedup.weight(3), 1u);

  dedup.addTrack(std::vector<Key>{0, 1, 3, 4, 5});

  // pixel doublet on the track: score 4 of 6
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1, 2}));
  // two pixels and one strip on the track: score 5 of 6
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1, 3, 9}));
  // strip seed, three of four strips on the track: score 3 of 4
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{3, 4, 5, 9}));
  // strip triplet, five of six strips on one track
  dedup.addTrack(std::vector<Key>{3, 4, 5, 6, 7});
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{3, 4, 5, 6, 7, 9}));

  // a seed of zero weight is never a duplicate
  dedup.setWeight(8, 0);
  dedup.addTrack(std::vector<Key>{8});
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{8}));
}

BOOST_AUTO_TEST_CASE(SeedDeduplicatorReset) {
  SeedDeduplicator dedup({});
  dedup.reset(5);
  dedup.setWeight(0, 3);
  dedup.addTrack(std::vector<Key>{0, 1});
  BOOST_CHECK(dedup.isDuplicate(std::vector<Key>{0, 1}));

  dedup.reset(5);
  BOOST_CHECK_EQUAL(dedup.nTracks(), 0u);
  BOOST_CHECK_EQUAL(dedup.weight(0), 1u);
  BOOST_CHECK(!dedup.isDuplicate(std::vector<Key>{0, 1}));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
