// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometryVisitor.hpp"

namespace Acts {

TrackingGeometryVisitor::~TrackingGeometryVisitor() = default;
TrackingGeometryMutableVisitor::~TrackingGeometryMutableVisitor() = default;

void TrackingGeometryVisitor::visitVolume(const TrackingVolume& /*volume*/) {}
void TrackingGeometryVisitor::visitPortal(const Portal& /*portal*/) {}
void TrackingGeometryVisitor::visitSurface(const Surface& /*surface*/) {}
void TrackingGeometryVisitor::visitLayer(const Layer& /*layer*/) {}
void TrackingGeometryVisitor::visitBoundarySurface(
    const BoundarySurfaceT<TrackingVolume>& /*boundary*/) {}

void TrackingGeometryMutableVisitor::visitVolume(TrackingVolume& /*volume*/) {}
void TrackingGeometryMutableVisitor::visitPortal(Portal& /*portal*/) {}
void TrackingGeometryMutableVisitor::visitSurface(Surface& /*surface*/) {}
void TrackingGeometryMutableVisitor::visitLayer(Layer& /*layer*/) {}
void TrackingGeometryMutableVisitor::visitBoundarySurface(
    BoundarySurfaceT<TrackingVolume>& /*boundary*/) {}

}  // namespace Acts
