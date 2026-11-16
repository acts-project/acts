// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

struct LinearizedTrack;
class Surface;
class MagneticFieldQuery;

/// Linearizer using a caller-owned field query. The query must not be retained.
/// Propagation uses the query's context and its own independent field cache.
using TrackLinearizer = Delegate<Result<LinearizedTrack>(
    const BoundTrackParameters& params, double linPointTime,
    const Surface& perigeeSurface, const GeometryContext& gctx,
    MagneticFieldQuery& field)>;

}  // namespace Acts
