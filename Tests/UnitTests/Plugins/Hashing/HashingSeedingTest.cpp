// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "ActsPlugins/Hashing/HashingAlgorithm.hpp"
#include "ActsPlugins/Hashing/HashingTraining.hpp"

#include <cstdlib>

#include <annoy/annoylib.h>
#include <annoy/kissrandom.h>

using namespace ActsPlugins;

namespace ActsTests {

// Function to create and initialize the test vector
Acts::SpacePointContainer2 createTestVector() {
  Acts::SpacePointContainer2 testVector(Acts::SpacePointColumns::X |
                                        Acts::SpacePointColumns::Y |
                                        Acts::SpacePointColumns::Z);
  testVector.reserve(6);

  const auto createSpacePoint = [&testVector](const float x, const float y,
                                              const float z) {
    auto sp = testVector.createSpacePoint();
    sp.x() = x;
    sp.y() = y;
    sp.z() = z;
    return sp;
  };

  createSpacePoint(27.2535, -18.0088, -146.526);
  createSpacePoint(42.9126, -27.3057, -171.477);
  createSpacePoint(42.7087, -27.4589, -171.557);
  createSpacePoint(74.3652, -45.8552, -221.277);
  createSpacePoint(104.12, -63.4203, -268.468);
  createSpacePoint(104.412, -63.1851, -268.468);

  return testVector;
}

BOOST_AUTO_TEST_SUITE(HashingSuite)

BOOST_AUTO_TEST_CASE(HashingBucketCreationTest) {
  // Initialize testVector using the createTestVector function
  auto testVector = createTestVector();

  /// Random seed for Annoy
  unsigned int annoySeed = 123456789;

  /// Number of features to use
  std::int32_t nf = 1;

  /// Size of the buckets = number of spacepoints in the bucket
  unsigned int bucketSize = 100;
  /// Number of zBins
  unsigned int zBins = 10000;
  /// Number of phiBins
  unsigned int phiBins = 0;

  /// Layer selection
  double layerRMin = 25;
  double layerRMax = 40;
  double layerZMin = -550;
  double layerZMax = 550;

  HashingAlgorithm::Config hashingConfig;
  hashingConfig.bucketSize = bucketSize;
  hashingConfig.zBins = zBins;
  hashingConfig.phiBins = phiBins;
  hashingConfig.layerRMin = layerRMin;
  hashingConfig.layerRMax = layerRMax;
  hashingConfig.layerZMin = layerZMin;
  hashingConfig.layerZMax = layerZMax;

  HashingTraining::Config hashingTrainingConfig;
  hashingTrainingConfig.annoySeed = annoySeed;
  hashingTrainingConfig.f = nf;

  HashingTraining hashingTraining(hashingTrainingConfig);
  HashingAlgorithm hashing(hashingConfig);

  // Hashing Training
  AnnoyModel annoyModel = hashingTraining.execute(testVector);

  // Hashing
  auto result = hashing.execute(annoyModel, testVector);

  // Check the number of buckets created
  BOOST_CHECK_GT(result.size(), 0);
}

BOOST_AUTO_TEST_CASE(HashingBucketContentTest) {
  // Initialize testVector using the createTestVector function
  auto testVector = createTestVector();

  /// Random seed for Annoy
  unsigned int annoySeed = 123456789;

  /// Number of features to use
  std::int32_t nf = 1;

  /// Size of the buckets = number of spacepoints in the bucket
  unsigned int bucketSize = 100;
  /// Number of zBins
  unsigned int zBins = 10000;
  /// Number of phiBins
  unsigned int phiBins = 0;

  /// Layer selection
  double layerRMin = 25;
  double layerRMax = 40;
  double layerZMin = -550;
  double layerZMax = 550;

  HashingAlgorithm::Config hashingConfig;
  hashingConfig.bucketSize = bucketSize;
  hashingConfig.zBins = zBins;
  hashingConfig.phiBins = phiBins;
  hashingConfig.layerRMin = layerRMin;
  hashingConfig.layerRMax = layerRMax;
  hashingConfig.layerZMin = layerZMin;
  hashingConfig.layerZMax = layerZMax;

  HashingTraining::Config hashingTrainingConfig;
  hashingTrainingConfig.annoySeed = annoySeed;
  hashingTrainingConfig.f = nf;

  HashingTraining hashingTraining(hashingTrainingConfig);
  HashingAlgorithm hashing(hashingConfig);

  // Hashing Training
  AnnoyModel annoyModel = hashingTraining.execute(testVector);

  // Hashing
  auto result = hashing.execute(annoyModel, testVector);

  // Validate bucket content
  for (const auto& bucket : result) {
    BOOST_CHECK_LE(bucket.size(), bucketSize);
    for (const auto& sp : bucket) {
      BOOST_CHECK_GE(sp, 0);
      BOOST_CHECK_LT(sp, testVector.size());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
