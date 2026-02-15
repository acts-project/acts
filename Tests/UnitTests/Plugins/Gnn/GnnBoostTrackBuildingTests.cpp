// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Helpers.hpp"
#include "ActsPlugins/Gnn/BoostTrackBuilding.hpp"

#include <algorithm>
#include <numeric>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GnnSuite)

BOOST_AUTO_TEST_CASE(test_track_building) {
  // Make some spacepoint IDs
  // The spacepoint ids are [100, 101, 102, ...]
  // They should not be zero based to check if the thing also works if the
  // spacepoint IDs do not match the node IDs used for the edges
  std::vector<int> spacePointIds(16);
  std::iota(spacePointIds.begin(), spacePointIds.end(), 100);

  // Build 4 tracks with 4 hits
  std::vector<std::vector<int>> refTracks;
  for (auto t = 0ul; t < 4; ++t) {
    refTracks.emplace_back(spacePointIds.begin() + 4 * t,
                           spacePointIds.begin() + 4 * (t + 1));
  }

  // Make edges
  std::vector<std::int64_t> e0, e1;
  for (const auto &track : refTracks) {
    for (auto it = track.begin(); it != track.end() - 1; ++it) {
      // edges must be 0 based, so subtract 100 again
      e0.push_back(*it - 100);
      e1.push_back(*std::next(it) - 100);
    }
  }

  ExecutionContext execCtx{Device::Cpu(), {}};

  auto edgeTensor = Tensor<std::int64_t>::Create({2, e0.size()}, execCtx);
  std::copy(e0.begin(), e0.end(), edgeTensor.data());
  std::copy(e1.begin(), e1.end(), edgeTensor.data() + e0.size());

  auto dummyNodes = Tensor<float>::Create({spacePointIds.size(), 16}, execCtx);
  auto dummyWeights = Tensor<float>::Create({e0.size(), 1}, execCtx);
  std::fill(dummyWeights.data(), dummyWeights.data() + dummyWeights.size(),
            1.f);

  // Run Track building
  auto logger = getDefaultLogger("TestLogger", Logging::ERROR);
  BoostTrackBuilding trackBuilder({}, std::move(logger));

  auto testTracks = trackBuilder({std::move(dummyNodes), std::move(edgeTensor),
                                  std::nullopt, std::move(dummyWeights)},
                                 spacePointIds);

  BOOST_CHECK_EQUAL(testTracks.size(), refTracks.size());

  // Sort tracks, so we can find them
  std::ranges::for_each(testTracks, [](auto &t) { std::ranges::sort(t); });
  std::ranges::sort(testTracks, std::less{}, [](auto &t) { return t.at(0); });

  std::ranges::for_each(refTracks, [](auto &t) { std::ranges::sort(t); });
  std::ranges::sort(refTracks, std::less{}, [](auto &t) { return t.at(0); });

  // Check what we have here
  for (const auto &refTrack : refTracks) {
    BOOST_CHECK(rangeContainsValue(testTracks, refTrack));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
