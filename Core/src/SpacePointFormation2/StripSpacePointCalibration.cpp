// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation2/StripSpacePointCalibration.hpp"

#include "Acts/SpacePointFormation2/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/detail/StdArrayLinalg.hpp"

Acts::OuterStripSpacePointCalibrationDetailsDerived
Acts::deriveOuterStripSpacePointCalibrationDetails(
    const OuterStripSpacePointCalibrationDetails& sp) {
  return detail::deriveOuterStripSpacePointCalibrationDetails(sp);
}

Eigen::Vector3f Acts::calibrateOuterStripSpacePoint(
    const Eigen::Vector3f& direction,
    const OuterStripSpacePointCalibrationDetailsDerived& sp) {
  std::array<float, 3> calibrated{};
  detail::calibrateOuterStripSpacePoint(detail::stdArrayCopy(direction), sp,
                                        calibrated,
                                        std::numeric_limits<float>::infinity());
  return detail::stdArrayToEigen(calibrated);
}
