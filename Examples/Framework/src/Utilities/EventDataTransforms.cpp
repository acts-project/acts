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

const ActsExamples::SimSpacePoint* ActsExamples::findSpacePointForIndex(
    ActsExamples::Index index, const SimSpacePointContainer& spacepoints) {
  auto match = [&](const SimSpacePoint& sp) {
    const auto& sls = sp.sourceLinks();
    return std::any_of(sls.begin(), sls.end(), [&](const auto& sl) {
      return sl.template get<IndexSourceLink>().index() == index;
    });
  };

  auto found = std::find_if(spacepoints.begin(), spacepoints.end(), match);

  if (found == spacepoints.end()) {
    return nullptr;
  }

  return &(*found);
}

ActsExamples::SimSeed ActsExamples::prototrackToSeed(
    const ActsExamples::ProtoTrack& track,
    const ActsExamples::SimSpacePointContainer& spacepoints) {
  auto findSpacePoint = [&](ActsExamples::Index index) {
    auto found = findSpacePointForIndex(index, spacepoints);
    if (found == nullptr) {
      throw std::runtime_error("No spacepoint found for source-link index " +
                               std::to_string(index));
    }
    return found;
  };

  const auto s = track.size();
  if (s < 3) {
    throw std::runtime_error(
        "Cannot convert track with less then 3 spacepoints to seed");
  }

  std::vector<const SimSpacePoint*> ps;
  ps.reserve(track.size());

  std::transform(track.begin(), track.end(), std::back_inserter(ps),
                 findSpacePoint);
  std::sort(ps.begin(), ps.end(),
            [](const auto& a, const auto& b) { return a->r() < b->r(); });

  // Simply use r = m*z + t and solve for r=0 to find z vertex position...
  // Probably not the textbook way to do
  const auto m =
      (ps.back()->r() - ps.front()->r()) / (ps.back()->z() - ps.front()->z());
  const auto t = ps.front()->r() - m * ps.front()->z();
  const auto z_vertex = -t / m;

  return SimSeed(*ps[0], *ps[s / 2], *ps[s - 1], z_vertex);
}
