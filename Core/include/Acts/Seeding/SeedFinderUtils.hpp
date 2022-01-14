// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"

namespace Acts {
/// @brief A partial description of a circle in u-v space.
struct LinCircle {
  float Zo;
  float cotTheta;
  float iDeltaR;
  float Er;
  float U;
  float V;
};

/// @brief Transform two spacepoints to a u-v space circle.
///
/// This function is a non-vectorized version of @a transformCoordinates.
///
/// @tparam external_spacepoint_t The external spacepoint type.
///
/// @param[in] sp The first spacepoint to use, either a bottom or top.
/// @param[in] spM The middle spacepoint to use.
/// @param[in] bottom Should be true if sp is a bottom SP.
template <typename external_spacepoint_t>
LinCircle transformCoordinates(
    const InternalSpacePoint<external_spacepoint_t>& sp,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom);

/// @brief Transform a vector of spacepoints to u-v space circles with respect
/// to a given middle spacepoint.
///
/// @tparam external_spacepoint_t The external spacepoint type.
///
/// @param[in] vec The list of bottom or top spacepoints
/// @param[in] spM The middle spacepoint.
/// @param[in] bottom Should be true if vec are bottom spacepoints.
/// @param[in] enableCutsForSortedSP enables sorting of cotTheta.
/// @param[out] linCircleVec The output vector to write to.
template <typename external_spacepoint_t>
void transformCoordinates(
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    bool enableCutsForSortedSP, std::vector<LinCircle>& linCircleVec);
}  // namespace Acts

#include "Acts/Seeding/SeedFinderUtils.ipp"
