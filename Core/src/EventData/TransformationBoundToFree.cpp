// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <algorithm>

Acts::FreeVector Acts::detail::transformBoundToFreeParameters(
    const Acts::Surface& surface, const GeometryContext& geoCtx,
    const Acts::BoundVector& boundParams) {
  // convert angles to global unit direction vector
  Vector3 direction = makeDirectionFromPhiTheta(boundParams[eBoundPhi],
                                                boundParams[eBoundTheta]);

  // convert local position to global position vector
  Vector2 local(boundParams[eBoundLoc0], boundParams[eBoundLoc1]);
  Vector3 position = surface.localToGlobal(geoCtx, local, direction);

  // construct full free-vector. time and q/p stay as-is.
  FreeVector freeParams = FreeVector::Zero();
  freeParams[eFreePos0] = position[ePos0];
  freeParams[eFreePos1] = position[ePos1];
  freeParams[eFreePos2] = position[ePos2];
  freeParams[eFreeTime] = boundParams[eBoundTime];
  freeParams[eFreeDir0] = direction[eMom0];
  freeParams[eFreeDir1] = direction[eMom1];
  freeParams[eFreeDir2] = direction[eMom2];
  freeParams[eFreeQOverP] = boundParams[eBoundQOverP];
  return freeParams;
}
