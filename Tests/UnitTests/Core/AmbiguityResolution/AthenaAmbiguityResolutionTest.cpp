// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/AmbiguityResolution/AthenaAmbiguityResolution.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/EventData/detail/TestTrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFinding/TrackSelector.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <map>

using Acts::MultiTrajectoryTraits::IndexType;

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(AthenaAmbiguityResolutionTest)

// Test fixture for AthenaAmbiguityResolution
struct Fixture {
  AthenaAmbiguityResolution::Config config;

  using TrackFeatures = AthenaAmbiguityResolution::TrackFeatures;
  using DetectorConfig = AthenaAmbiguityResolution::DetectorConfig;

  Fixture() {
    // Set up any resources used by the tests
    config.volumeMap = {{8, 0},  {9, 0},  {10, 0}, {13, 0}, {14, 0},
                        {15, 0}, {16, 0}, {17, 0}, {18, 0}, {19, 0},
                        {20, 0}, {22, 1}, {23, 1}, {24, 1}};

    auto tempDetector = DetectorConfig();
    std::map<std::size_t, DetectorConfig> detectorMap = {{0, tempDetector},
                                                         {1, tempDetector}};
    config.detectorMap = detectorMap;

    config.minScore = 0;
    config.minScoreSharedTracks = 100;
    config.maxShared = 5;
    config.maxSharedTracksPerMeasurement = 10;
    config.phiMax = 3.14;
    config.phiMin = -3.14;
    config.etaMax = 2.7;
    config.etaMin = -2.7;
    config.pTMax = 1400;
    config.pTMin = 0.5;
    config.useAmbiguityFunction = false;
  }

  ~Fixture() = default;
};

// Helper function to create a sample input for getCleanedOutTracks
std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
createSampleInput() {
  std::vector<std::pair<std::size_t, std::vector<std::size_t>>> trackVolumes = {
      {0, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
      {1, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10}},
      {2, {13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8}},
      {3, {13, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8}},
      {4, {19, 18, 18, 18, 10, 10, 10, 10, 10, 10, 10, 10, 10}}};

  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack;
  // Add sample measurements for each track

  for (const auto& trackVolume : trackVolumes) {
    std::vector<std::tuple<std::size_t, std::size_t, bool>> measurements;
    for (std::size_t i = 0; i < trackVolume.second.size(); ++i) {
      measurements.push_back(
          std::make_tuple(i + 2, trackVolume.second[i], false));
    }
    measurementsPerTrack.push_back(measurements);
  }

  return measurementsPerTrack;
}

BOOST_FIXTURE_TEST_CASE(GetCleanedOutTracksTest, Fixture) {
  Fixture fixture;
  // Create an instance of AthenaAmbiguityResolution
  AthenaAmbiguityResolution tester(fixture.config);

  // Create sample input
  std::vector<std::vector<std::tuple<std::size_t, std::size_t, bool>>>
      measurementsPerTrack = createSampleInput();

  std::vector<double> TrackSore;
  for (std::size_t i = 0; i < measurementsPerTrack.size(); i++) {
    TrackSore.push_back(60 + 40 * i);
  }

  std::vector<std::map<std::size_t, TrackFeatures>> trackFeaturesMaps = {
      {{0, {0, 14, 0, 0}}, {1, {0, 0, 0, 0}}},
      {{0, {0, 15, 0, 0}}, {1, {0, 0, 0, 0}}},
      {{0, {0, 17, 0, 0}}, {1, {0, 0, 0, 0}}},
      {{0, {0, 18, 0, 0}}, {1, {0, 0, 0, 0}}},
      {{0, {0, 14, 0, 0}}, {1, {0, 0, 0, 0}}}};

  // Call the function under testBOOST_FIXTURE_TEST_CASE
  std::vector<std::size_t> cleanTracks = tester.getCleanedOutTracks(
      TrackSore, trackFeaturesMaps, measurementsPerTrack);

  // Assert the expected results
  BOOST_CHECK_EQUAL(measurementsPerTrack.size(), 5);
  BOOST_CHECK_EQUAL(cleanTracks.size(), 3);

  for (std::size_t i = 0; i < cleanTracks.size(); i++) {
    auto score = TrackSore[cleanTracks[i]];
    BOOST_CHECK_GT(score, fixture.config.minScoreSharedTracks);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
