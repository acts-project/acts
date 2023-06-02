// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/EventDataTransforms.hpp"

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

ActsExamples::SimSeed ActsExamples::ProtoTrackToSeed::operator()(
    const ActsExamples::ProtoTrack& track) const {
  const auto s = track.size();
  if (s < 3) {
    throw std::runtime_error(
        "Cannot convert track with less then 3 spacepoints to seed");
  }

  std::array<SimSpacePoint, 3> ps{{m_spacePoints.at(track[0]),
                                   m_spacePoints.at(track[s / 2]),
                                   m_spacePoints.at(track[s - 1])}};
  std::sort(ps.begin(), ps.end(),
            [](const auto& a, const auto& b) { return a.r() < b.r(); });

  // Simply use r = m*z + t and solve for r=0 to find z vertex position...
  // Probably not the textbook way to do
  const auto m =
      (ps.back().r() - ps.front().r()) / (ps.back().z() - ps.front().z());
  const auto t = ps.front().r() - m * ps.front().z();
  const auto z_vertex = -t / m;

  return SimSeed(ps[0], ps[1], ps[2], z_vertex);
}
