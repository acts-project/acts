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
#include "Acts/Utilities/Logger.hpp"

Acts::BoundVector Acts::detail::transformFreeToBoundParameters(
    const FreeVector& freeParams, const Surface& surface,
    const GeometryContext& geoCtx) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto position = freeParams.segment<3>(eFreePos0);
  auto direction = freeParams.segment<3>(eFreeDir0);
  auto result = surface.globalToLocal(geoCtx, position, direction);
  if (result.ok()) {
    auto localPosition = result.value();
    bp[eBoundLoc0] = localPosition[ePos0];
    bp[eBoundLoc1] = localPosition[ePos1];
  } else {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("ParameterTransformation", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation from free to bound.")
  }
  bp[eBoundTime] = freeParams[eFreeTime];
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = freeParams[eFreeQOverP];
  return bp;
}

Acts::BoundVector Acts::detail::transformFreeToBoundParameters(
    const Acts::Vector3D& position, FreeScalar time,
    const Acts::Vector3D& direction, FreeScalar qOverP,
    const Acts::Surface& surface, const Acts::GeometryContext& geoCtx) {
  // initialize the bound vector
  BoundVector bp = BoundVector::Zero();
  // convert global to local position on the surface
  auto result = surface.globalToLocal(geoCtx, position, direction);
  if (result.ok()) {
    auto localPosition = result.value();
    bp[eBoundLoc0] = localPosition[ePos0];
    bp[eBoundLoc1] = localPosition[ePos1];
  } else {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("ParameterTransformation", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation from free to bound.")
  }
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}

Acts::BoundVector Acts::detail::transformFreeToCurvilinearParameters(
    FreeScalar time, FreeScalar phi, FreeScalar theta, FreeScalar qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = phi;
  bp[eBoundTheta] = theta;
  bp[eBoundQOverP] = qOverP;
  return bp;
}

Acts::BoundVector Acts::detail::transformFreeToCurvilinearParameters(
    FreeScalar time, const Vector3D& direction, FreeScalar qOverP) {
  BoundVector bp = BoundVector::Zero();
  // local coordinates are zero by construction
  bp[eBoundTime] = time;
  bp[eBoundPhi] = VectorHelpers::phi(direction);
  bp[eBoundTheta] = VectorHelpers::theta(direction);
  bp[eBoundQOverP] = qOverP;
  return bp;
}
