// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

namespace Acts::Experimental::detail {

template <CompositeSpacePoint Point_t>
int CompSpacePointAuxiliaries::strawSign(const Line_t& line,
                                         const Point_t& strawSp) {
  return strawSign(line.position(), line.direction(), strawSp);
}

template <CompositeSpacePoint Point_t>
int CompSpacePointAuxiliaries::strawSign(const Vector& pos, const Vector& dir,
                                         const Point_t& strawSp) {
  if (!strawSp.isStraw()) {
    return 0;
  }
  const double dist = Acts::detail::LineHelper::signedDistance(
      pos, dir, , strawSp.localPosition(), strawSp.sensorDirection());
  return dist > 0. ? 1 : -1;
}
template <CompositeSpacePointContainer StrawCont_t>
std::vector<int> CompSpacePointAuxiliaries::strawSigns(
    const Line_t& line, const StrawCont_t& measurements) {
  return strawSigns(line.position(), line.direction(), measurements);
}
template <CompositeSpacePointContainer StrawCont_t>
std::vector<int> CompSpacePointAuxiliaries::strawSigns(
    const Vector& pos, const Vector& dir, const StrawCont_t& measurements) {
  signs.reserve(measurements.size());
  for (const auto& strawSp : measurements) {
    signs.push_back(strawSign(pos, dir, *strawSp));
  }
  return signs;
}
template <CompositeSpacePoint Point_t>
void CompSpacePointAuxiliaries::updateSpatialResidual(
    const Line_t& line, const Point_t& spacePoint) {
  if (spacePoint.isStraw()) {
    /// Fetch the hit position & direction
    const auto& wireDir{spacePoint.sensorDirection()};
    /// Calculate the distance from the two reference points
    const Vector hitMinSeg = spacePoint.localPosition() - line.position();

    if (!updateStrawResidual(line, hitMinSeg, wireDir,
                             spacePoint.driftRadius())) {
      return;
    }

    if (m_cfg.calcAlongStraw && spacePoint.measuresLoc0()) {
      /// If the tube is a twin-tube, the hit position is no longer arbitrary
      /// along the wire. Calculate the distance along the wire towards the
      /// point of closest approach.
      updateAlongTheStraw(line, hitMinSeg, wireDir);
    }
  } else {
    updateStripResidual(line, spacePoint.planeNormal(),
                        spacePoint.toNextSensor(), spacePoint.sensorDirection(),
                        spacePoint.localPosition(), spacePoint.measuresLoc1(),
                        spacePoint.measuresLoc0());
  }
}

template <CompositeSpacePoint Point_t>
void CompSpacePointAuxiliaries::updateFullResidual(const Line_t& line,
                                                   const double timeOffset,
                                                   const Point_t& spacePoint,
                                                   const double driftV,
                                                   const double driftA) {
  /// Calculate first the spatial residual
  updateSpatialResidual(line, spacePoint);

  /// Calculate the time residual for strip-like measurements
  if (!spacePoint.isStraw()) {
    /// If the measurement does not provide time, then simply reset the time
    /// partial components
    if (!spacePoint.hasTime()) {
      resetTime();
      return;
    }
    updateTimeStripRes(spacePoint.toNextSensor(), spacePoint.sensorDirection(),
                       spacePoint.localPosition(), spacePoint.measuresLoc1(),
                       spacePoint.time(), timeOffset);
  } else {
    updateTimeStrawRes(line, spacePoint.localPosition() - line.position(),
                       spacePoint.sensorDirection(), spacePoint.driftRadius(),
                       driftV, driftA);
  }
}

}  // namespace Acts::Experimental::detail
