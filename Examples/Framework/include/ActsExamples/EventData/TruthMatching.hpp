// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <cstdint>
#include <map>
#include <optional>

namespace ActsExamples {

struct TrackMatchEntry {
  std::optional<SimBarcode> particle;
};

struct ParticleMatchEntry {
  std::optional<TrackIndex> track;
  std::uint32_t duplicates{};
  std::uint32_t fakes{};
};

using TrackParticleMatching = std::map<TrackIndex, TrackMatchEntry>;
using ParticleTrackMatching = std::map<SimBarcode, ParticleMatchEntry>;

}  // namespace ActsExamples
