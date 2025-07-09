// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <cstdint>
#include <limits>

namespace Acts {

using SpacePointIndex2 = std::uint32_t;
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;

using SeedIndex2 = std::uint32_t;

using TrackIndexType = std::uint32_t;
static constexpr TrackIndexType kTrackIndexInvalid =
    std::numeric_limits<TrackIndexType>::max();

using ProjectorBitset = std::uint64_t;

template <std::size_t measdim>
using SubspaceIndices = std::array<std::uint8_t, measdim>;
using BoundSubspaceIndices = SubspaceIndices<eBoundSize>;
static constexpr BoundSubspaceIndices kBoundSubspaceIndicesInvalid = {
    eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize};
using SerializedSubspaceIndices = std::uint64_t;

}  // namespace Acts
