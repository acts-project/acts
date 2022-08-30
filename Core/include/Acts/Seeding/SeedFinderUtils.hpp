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
#include "Acts/Seeding/SeedfinderConfig.hpp"

namespace Acts {
/// @brief A partial description of a circle in u-v space.
struct LinCircle {
  float Zo;
  float cotTheta;
  float iDeltaR;
  float Er;
  float U;
  float V;
  float x;
  float y;
  float z;
  float r;
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
LinCircle transformCoordinates(InternalSpacePoint<external_spacepoint_t>& sp,
                               InternalSpacePoint<external_spacepoint_t>& spM,
                               bool bottom);

template <typename external_spacepoint_t, typename callable_t>
LinCircle transformCoordinates(external_spacepoint_t& sp,
                               external_spacepoint_t& spM, bool bottom,
                               callable_t&& extractFunction);

/// @brief Transform a vector of spacepoints to u-v space circles with respect
/// to a given middle spacepoint.
///
/// @tparam external_spacepoint_t The external spacepoint type.
///
/// @param[in] vec The list of bottom or top spacepoints
/// @param[in] spM The middle spacepoint.
/// @param[in] bottom Should be true if vec are bottom spacepoints.
/// @param[out] linCircleVec The output vector to write to.
template <typename external_spacepoint_t>
void transformCoordinates(
    std::vector<InternalSpacePoint<external_spacepoint_t>*>& vec,
    InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec);

template <typename external_spacepoint_t, typename callable_t>
void transformCoordinates(std::vector<external_spacepoint_t*>& vec,
                          external_spacepoint_t& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec,
                          callable_t&& extractFunction);

/// @brief Check the compatibility of spacepoint coordinates in xyz assuming the Bottom-Middle direction with the strip meassument details
///
/// @tparam external_spacepoint_t The external spacepoint type.
/// @tparam sp_range_t Container type for the space point collections.
///
/// @param[in] config Seedfinder config containing the delegates to the strip measurement details.
/// @param[in] sp Input space point used in the check.
/// @param[in] spacepointPosition Spacepoint coordinates in xyz plane.
/// @param[in] toleranceParam Parameter used to evaluate if spacepointPosition is inside the detector elements.
/// @param[out] outputCoordinates The output vector to write to.
/// @returns Boolean that says if spacepoint is compatible with being inside the detector element.
template <typename external_spacepoint_t, typename sp_range_t>
bool xyzCoordinateCheck(Acts::SeedfinderConfig<external_spacepoint_t> config,
                        sp_range_t sp, const double* spacepointPosition,
                        const float toleranceParam, double* outputCoordinates);

}  // namespace Acts

#include "Acts/Seeding/SeedFinderUtils.ipp"
