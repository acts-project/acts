// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"

class GeoFullPhysVol;

namespace Acts {
/// Collect the sensitive surface & detector element
using GeoModelSensitiveSurface =
    std::tuple<std::shared_ptr<GeoModelDetectorElement>,
               std::shared_ptr<Surface>>;

namespace GeoModelSurfaceConverter {
/// @brief conversion to sensitive surface
///
/// @param geoPhysVol the geoPhysVol to convert
///
/// @return a detector element and a surface
GeoModelSensitiveSurface convertToSensitiveSurface(
    const GeoFullPhysVol& geoPhysVol);
}  // namespace GeoModelSurfaceConverter
}  // namespace Acts
