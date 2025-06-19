// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/ExaTrkX/BoostTrackBuilding.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <numeric>

BOOST_AUTO_TEST_CASE(test_track_building) {
  // Make some spacepoint IDs
  // The spacepoint ids are [100, 101, 102, ...]
  // They should not be zero based to check if the thing also works if the
  // spacepoint IDs do not match the node IDs used for the edges
  std::vector<int> spacepointIds(16);
  std::iota(spacepointIds.begin(), spacepointIds.end(), 100);

  // Build 4 tracks with 4 hits
  std::vector<std::vector<int>> refTracks;
  for (auto t = 0ul; t < 4; ++t) {
    refTracks.emplace_back(spacepointIds.begin() + 4 * t,
                           spacepointIds.begin() + 4 * (t + 1));
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

  Acts::ExecutionContext execCtx{Acts::Device::Cpu(), {}};

  auto edgeTensor = Acts::Tensor<std::int64_t>::Create({2, e0.size()}, execCtx);
  std::copy(e0.begin(), e0.end(), edgeTensor.data());
  std::copy(e1.begin(), e1.end(), edgeTensor.data() + e0.size());

  auto dummyNodes =
      Acts::Tensor<float>::Create({spacepointIds.size(), 16}, execCtx);
  auto dummyWeights = Acts::Tensor<float>::Create({e0.size(), 1}, execCtx);
  std::fill(dummyWeights.data(), dummyWeights.data() + dummyWeights.size(),
            1.f);

  // Run Track building
  auto logger = Acts::getDefaultLogger("TestLogger", Acts::Logging::ERROR);
  Acts::BoostTrackBuilding trackBuilder({}, std::move(logger));

  auto testTracks = trackBuilder({std::move(dummyNodes), std::move(edgeTensor),
                                  std::nullopt, std::move(dummyWeights)},
                                 spacepointIds);

  BOOST_CHECK_EQUAL(testTracks.size(), refTracks.size());

  // Sort tracks, so we can find them
  std::ranges::for_each(testTracks, [](auto &t) { std::ranges::sort(t); });
  std::ranges::sort(testTracks, std::less{}, [](auto &t) { return t.at(0); });

  std::ranges::for_each(refTracks, [](auto &t) { std::ranges::sort(t); });
  std::ranges::sort(refTracks, std::less{}, [](auto &t) { return t.at(0); });

  // Check what we have here
  for (const auto &refTrack : refTracks) {
    BOOST_CHECK(Acts::rangeContainsValue(testTracks, refTrack));
  }
}
