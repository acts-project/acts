// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/CylinderSurface.hpp"

namespace ActsTests {
/// Prototype of a full-cylinder on-axis fast path. Retains both intersections
/// and the existing bounds/status semantics. All other cases use Core.
Acts::MultiIntersection3D intersectOnAxisCylinder(
    const Acts::CylinderSurface& surface, const Acts::GeometryContext& gctx,
    const Acts::Vector3& position, const Acts::Vector3& direction,
    const Acts::BoundaryTolerance& boundaryTolerance,
    double tolerance = Acts::s_onSurfaceTolerance);
}  // namespace ActsTests
