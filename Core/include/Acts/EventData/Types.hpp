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
#include <span>

namespace Acts {

/// Type alias for track index values
using TrackIndexType = std::uint32_t;
static constexpr TrackIndexType kTrackIndexInvalid =
    std::numeric_limits<TrackIndexType>::max();

/// @brief Type alias for bitset representing parameter projections
/// @details Used to store parameter projection information in a compact bit format
using ProjectorBitset = std::uint64_t;

/// Type alias for subspace index
using SubspaceIndex = std::uint8_t;
/// Template alias for subspace indices array
template <std::size_t measdim>
using SubspaceIndices = std::array<SubspaceIndex, measdim>;
/// Type alias for bound parameter subspace indices
using BoundSubspaceIndices = SubspaceIndices<eBoundSize>;

// template <std::size_t measdim>
/// @brief type alias for indices of bound track parameters subspace
/// @details used to specify which components of the bound track parameters are being referenced
// using BoundSubspaceIndices = std::array<std::size_t, eBoundSize>;
static constexpr BoundSubspaceIndices kBoundSubspaceIndicesInvalid = {
    eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize, eBoundSize};

/// @brief Type alias for serialized subspace indices
/// @details Compact representation of subspace indices as a 64-bit unsigned integer
using SerializedSubspaceIndices = std::uint64_t;

using SpacePointIndex2 = std::uint32_t;
using SpacePointIndexRange2 = std::pair<SpacePointIndex2, SpacePointIndex2>;
using SpacePointIndexSubset2 = std::span<const SpacePointIndex2>;

using SeedIndex2 = std::uint32_t;

}  // namespace Acts
