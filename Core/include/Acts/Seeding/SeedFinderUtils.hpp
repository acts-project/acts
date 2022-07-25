// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePoint.hpp"

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
/// @param[in] sp The first spacepoint to use, either a bottom or top.
/// @param[in] spM The middle spacepoint to use.
/// @param[in] bottom Should be true if sp is a bottom SP.
LinCircle transformCoordinates(Acts::SpacePoint& sp, Acts::SpacePoint& spM,
                               bool bottom);

/// @brief Transform a vector of spacepoints to u-v space circles with respect
/// to a given middle spacepoint.
///
/// @param[in] vec The list of bottom or top spacepoints
/// @param[in] spM The middle spacepoint.
/// @param[in] bottom Should be true if vec are bottom spacepoints.
/// @param[out] linCircleVec The output vector to write to.
void transformCoordinates(std::vector<Acts::SpacePoint*>& vec,
                          Acts::SpacePoint& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec);

void transformCoordinates(std::vector<Acts::SpacePoint*>& vec,
                          Acts::SpacePoint& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec);

/// @brief Check the compatibility of spacepoint coordinates in xyz assuming the Bottom-Middle direction with the strip measurement details
///
/// @tparam sp_range_t Container type for the space point collections.
///
/// @param[in] config Seedfinder config containing the delegates to the strip measurement details.
/// @param[in] sp Input space point used in the check.
/// @param[in] spacepointPosition Spacepoint coordinates in xyz plane.
/// @param[in] toleranceParam Parameter used to evaluate if spacepointPosition is inside the detector elements.
/// @param[out] outputCoordinates The output vector to write to.
/// @returns Boolean that says if spacepoint is compatible with being inside the detector element.
bool xyzCoordinateCheck(Acts::SeedfinderConfig config, Acts::SpacePoint* sp,
                        const double* spacepointPosition,
                        const float toleranceParam, double* outputCoordinates);

}  // namespace Acts