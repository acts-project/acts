// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include <vector>

ActsExamples::ProtoTrack ActsExamples::seedToPrototrack(
    const ActsExamples::SimSeed& seed) {
  ProtoTrack track;
  track.reserve(seed.sp().size());
  for (auto spacePointPtr : seed.sp()) {
    for (const auto& slink : spacePointPtr->sourceLinks()) {
      const auto& islink = slink.get<IndexSourceLink>();
      track.emplace_back(islink.index());
    }
  }
  return track;
}

ActsExamples::ProtoTrackContainer ActsExamples::seedsToPrototracks(
    const ActsExamples::SimSeedContainer& seeds) {
  ProtoTrackContainer tracks;
  tracks.reserve(seeds.size());

  for (const auto& seed : seeds) {
    tracks.push_back(seedToPrototrack(seed));
  }

  return tracks;
}
