// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"

#include "GeoModelKernel/GeoDefinitions.h"
class GeoShape;

namespace Acts::GeoModel {

/// @brief Convert a GeoModel shape to a DetectorVolume
///
/// @param shape the GeoModel shape
/// @param transform the transform to be applied
/// @return the DetectorVolume
std::shared_ptr<Experimental::DetectorVolume> convertVolume(
    const GeometryContext& context, const GeoShape& shape,
    const std::string& name, const GeoTrf::Transform3D& transform,
    const std::vector<GeoModelSensitiveSurface>& sensitives);

}  // namespace Acts::GeoModel
