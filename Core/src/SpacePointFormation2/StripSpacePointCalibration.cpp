// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation2/StripSpacePointCalibration.hpp"

#include "Acts/SpacePointFormation2/detail/StripSpacePointCalibrationImpl.hpp"

Acts::StripSpacePointCalibrationDetailsDerived
Acts::deriveStripSpacePointCalibrationDetails(
    const StripSpacePointCalibrationDetails& sp) {
  return detail::deriveStripSpacePointCalibrationDetails(sp);
}

std::optional<Eigen::Vector3f> Acts::calibrateStripSpacePoint(
    const StripSpacePointCalibrationDetailsDerived& sp,
    const Eigen::Vector3f& direction, float tolerance) {
  Eigen::Vector3f calibrated;
  if (!detail::calibrateStripSpacePoint(
          sp, std::span<const float, 3>(direction.data(), direction.size()),
          std::span<float, 3>(calibrated.data(), calibrated.size()),
          tolerance)) {
    return std::nullopt;
  }
  return calibrated;
}
