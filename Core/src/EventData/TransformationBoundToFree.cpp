// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

Acts::FreeVector Acts::detail::transformBoundToFreeParameters(
    const Acts::Surface& surface, const GeometryContext& geoCtx,
    const Acts::BoundVector& boundParams) {
  // convert angles to global unit direction vector
  const auto dir = makeDirectionUnitFromPhiTheta(boundParams[eBoundPhi],
                                                 boundParams[eBoundTheta]);

  // convert local position to global position vector
  Vector2D loc(boundParams[eBoundLoc0], boundParams[eBoundLoc1]);
  Vector3D pos = Vector3D::Zero();
  surface.localToGlobal(geoCtx, loc, dir, pos);

  // construct full free-vector. time and q/p stay as-is.
  FreeVector freeParams = FreeVector::Zero();
  freeParams[eFreePos0] = pos[ePos0];
  freeParams[eFreePos1] = pos[ePos1];
  freeParams[eFreePos2] = pos[ePos2];
  freeParams[eFreeTime] = boundParams[eBoundTime];
  freeParams[eFreeDir0] = dir[eMom0];
  freeParams[eFreeDir1] = dir[eMom1];
  freeParams[eFreeDir2] = dir[eMom2];
  freeParams[eFreeQOverP] = boundParams[eBoundQOverP];
  return freeParams;
}
