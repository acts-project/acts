// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"

#include <optional>

#include <Eigen/Core>

namespace Acts {

StripSpacePointCalibrationDetailsDerived
deriveStripSpacePointCalibrationDetails(
    const StripSpacePointCalibrationDetails& sp);

std::optional<Eigen::Vector3f> calibrateStripSpacePoint(
    const StripSpacePointCalibrationDetailsDerived& sp,
    const Eigen::Vector3f& direction, float tolerance);

}  // namespace Acts
