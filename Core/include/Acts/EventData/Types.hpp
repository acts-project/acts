// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <cstdint>
#include <limits>

namespace Acts {

using TrackIndexType = std::uint32_t;
static constexpr TrackIndexType kTrackIndexInvalid =
    std::numeric_limits<TrackIndexType>::max();

using ProjectorBitset = std::uint64_t;

template <std::size_t measdim>
using ProjectorMapping = std::array<std::uint8_t, measdim>;
using FullProjectorMapping = ProjectorMapping<eBoundSize>;
static constexpr FullProjectorMapping kFullProjectorMappingEmpty = {
    eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize};

}  // namespace Acts
