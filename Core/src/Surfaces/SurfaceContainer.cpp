// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceContainer.hpp"

Acts::SurfacePtrsContainer Acts::SurfaceContainer::getPtrs(
    const DetectorPtr& detector) const {
    Acts::SurfaceContainer::SurfaceVisitor visitor;
    detector->visitSurfaces(visitor);
    return visitor.surfacePtrs;
}

Acts::SurfacePtrsContainer Acts::SurfaceContainer::getPtrs(
    const TrackingGeometryPtr& tGeometryPtr) const {
  Acts::SurfaceContainer::SurfaceVisitor visitor;
  tGeometryPtr->visitSurfaces(visitor);
  return visitor.surfacePtrs;
}
