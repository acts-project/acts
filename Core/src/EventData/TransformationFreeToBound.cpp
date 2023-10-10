// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationFreeToBound.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>

Acts::Result<Acts::BoundVector> Acts::detail::transformFreeToBoundParameters(
    const FreeVector& freeParams, const Surface& surface,
    const GeometryContext& geoCtx, ActsScalar tolerance) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto position = freeParams.segment<3>(eFreePos0);
  auto direction = freeParams.segment<3>(eFreeDir0);
  auto result = surface.globalToLocal(geoCtx, position, direction, tolerance);
  if (!result.ok()) {
    return Result<Acts::BoundVector>::failure(result.error());
  }

  auto localPosition = result.value();
  bp[eBoundLoc0] = localPosition[ePos0];
  bp[eBoundLoc1] = localPosition[ePos1];

  bp[eBoundTime] = freeParams[eFreeTime];
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = freeParams[eFreeQOverP];
  return Result<Acts::BoundVector>::success(bp);
}

Acts::Result<Acts::BoundVector> Acts::detail::transformFreeToBoundParameters(
    const Acts::Vector3& position, ActsScalar time,
    const Acts::Vector3& direction, ActsScalar qOverP,
    const Acts::Surface& surface, const Acts::GeometryContext& geoCtx,
    ActsScalar tolerance) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto result = surface.globalToLocal(geoCtx, position, direction, tolerance);
  if (!result.ok()) {
    return Result<Acts::BoundVector>::failure(result.error());
  }

  auto localPosition = result.value();
  bp[eBoundLoc0] = localPosition[ePos0];
  bp[eBoundLoc1] = localPosition[ePos1];

  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return Result<Acts::BoundVector>::success(bp);
}

Acts::BoundVector Acts::detail::transformFreeToCurvilinearParameters(
    ActsScalar time, ActsScalar phi, ActsScalar theta, ActsScalar qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = phi;
  bp[eBoundTheta] = theta;
  bp[eBoundQOverP] = qOverP;
  return bp;
}

Acts::BoundVector Acts::detail::transformFreeToCurvilinearParameters(
    ActsScalar time, const Vector3& direction, ActsScalar qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}
