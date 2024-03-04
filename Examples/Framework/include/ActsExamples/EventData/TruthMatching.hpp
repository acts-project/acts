// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <cstdint>
#include <map>
#include <optional>

namespace ActsExamples {

struct TrackMatchEntry {
  std::optional<SimBarcode> particle;

  /// Number of hits on the track that are associated to a particle
  /// Sorted by decreasing number of hits
  std::vector<ParticleHitCount> contributingParticles;

  bool isDuplicate{};
  bool isFake{};
};

struct ParticleMatchEntry {
  std::optional<TrackIndexType> track;
  std::uint32_t duplicates{};
  std::uint32_t fakes{};
};

using ProtoTrackParticleMatching = std::map<TrackIndexType, TrackMatchEntry>;
using ParticleProtoTrackMatching = std::map<SimBarcode, ParticleMatchEntry>;

using TrackParticleMatching = std::map<TrackIndexType, TrackMatchEntry>;
using ParticleTrackMatching = std::map<SimBarcode, ParticleMatchEntry>;

}  // namespace ActsExamples
