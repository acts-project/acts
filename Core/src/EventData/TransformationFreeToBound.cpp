// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationFreeToBound.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::BoundVector Acts::detail::transformFreeToBoundParameters(
    const Acts::Vector3D& position, double time,
    const Acts::Vector3D& direction, double qOverP,
    const Acts::Surface& surface, const Acts::GeometryContext& geoCtx) {
  // this assumes the position is already on the surface
  Vector2D localPosition;
  surface.globalToLocal(geoCtx, position, direction, localPosition);
  // construct the bound vector
  BoundVector bp = BoundVector::Zero();
  bp[eBoundLoc0] = localPosition[ePos0];
  bp[eBoundLoc1] = localPosition[ePos1];
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}

Acts::BoundVector Acts::detail::transformFreeToCurvilinearParameters(
    double time, const Acts::Vector3D& direction, double qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}
