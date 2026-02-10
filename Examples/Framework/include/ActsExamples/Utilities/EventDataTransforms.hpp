// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"

namespace ActsExamples {

ProtoTrack seedToPrototrack(const ConstSeedProxy &seed);

std::optional<ConstSpacePointProxy> findSpacePointForIndex(
    Index index, const SpacePointContainer &spacePoints);

SeedProxy prototrackToSeed(const ProtoTrack &track,
                           const SpacePointContainer &spacePoints,
                           SeedContainer &seeds);

}  // namespace ActsExamples
