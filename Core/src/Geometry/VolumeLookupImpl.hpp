// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
class TrackingVolume;
namespace detail {
// Select direct children linearly and descend through their configured policy.
// Accelerator fallbacks bypass dispatch on the root to avoid recursion cycles.
const TrackingVolume* linearVolumeScan(const TrackingVolume& root,
                                       const GeometryContext& gctx,
                                       const Vector3& position,
                                       double tolerance);

}  // namespace detail
}  // namespace Acts
