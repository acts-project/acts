// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"

#include <Eigen/Core>

namespace Acts {

/// Derives the strip space point calibration details.
/// @param sp The strip space point calibration details
/// @return The derived strip space point calibration details
OuterStripSpacePointCalibrationDetailsDerived
deriveOuterStripSpacePointCalibrationDetails(
    const OuterStripSpacePointCalibrationDetails& sp);

/// Calibrates the strip space point using the assumed particle direction and
/// the strip space point calibration details.
/// @note This function does not check if the calibrated space point lies within the strip.
/// @param sp The strip space point calibration details
/// @param direction The assumed particle direction
/// @return The calibrated outer strip space point
Eigen::Vector3f calibrateOuterStripSpacePoint(
    const Eigen::Vector3f& direction,
    const OuterStripSpacePointCalibrationDetailsDerived& sp);

}  // namespace Acts
