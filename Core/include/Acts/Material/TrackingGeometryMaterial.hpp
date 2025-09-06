// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"

#include <map>
#include <memory>
#include <utility>

namespace Acts {

/// Type alias for surface material maps indexed by geometry identifier
using SurfaceMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
/// Type alias for volume material maps indexed by geometry identifier
using VolumeMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
/// Type alias for tracking geometry material containing surface and volume maps
using TrackingGeometryMaterial =
    std::pair<SurfaceMaterialMaps, VolumeMaterialMaps>;

}  // namespace Acts
