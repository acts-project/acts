// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

namespace Acts::detail {

template <StationSpacePoint Point_t>
int StrawLineFitAuxiliaries::strawSign(const Line_t& line,
                                       const Point_t& strawSp) {
  if (!strawSp.isStraw()) {
    return 0;
  }
  const double dist = LineHelper::signedDistance(
      line.position(), line.direction(), strawSp.localPosition(),
      strawSp.sensorDirection());
  return dist > 0. ? 1 : -1;
}

template <StationSpacePoint Point_t>
void StrawLineFitAuxiliaries::updateSpatialResidual(const Line_t& line,
                                                    const Point_t& spacePoint) {
  if (spacePoint.isStraw()) {
    /// Fetch the hit position & direction
    const auto& wireDir{spacePoint.sensorDirection()};
    /// Calculate the distance from the two reference points
    const Vector hitMinSeg = spacePoint.localPosition() - line.position();

    if (!updateStrawResidual(line, hitMinSeg, wireDir,
                             spacePoint.driftRadius())) {
      return;
    }

    if (m_cfg.calcAlongStraw && spacePoint.measNonPrecCoord()) {
      /// If the tube is a twin-tube, the hit position is no longer arbitrary
      /// along the wire. Calculate the distance along the wire towards the
      /// point of closest approach.
      updateAlongTheStraw(line, hitMinSeg, wireDir);
    }
  } else {
    updateStripResidual(line, spacePoint.planeNormal(),
                        spacePoint.sensorNormal(), spacePoint.sensorDirection(),
                        spacePoint.localPosition(), spacePoint.measPrecCoord(),
                        spacePoint.measNonPrecCoord());
  }
}

}  // namespace Acts::detail
