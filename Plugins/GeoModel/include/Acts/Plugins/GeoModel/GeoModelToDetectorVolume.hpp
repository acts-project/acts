// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"

#include "GeoModelKernel/GeoDefinitions.h"
class GeoShape;

namespace Acts::GeoModel {

// function converts GeoShape to Acts::Volume
Volume convertVolume(const Transform3& trf, const GeoShape& shape);

/// @brief Convert a GeoModel shape to a DetectorVolume
///
/// @param shape the GeoModel shape
/// @param transform the transform to be applied
/// @param name the name of the volume
/// @param sensitives the sensitive surfaces
/// @return the pair of Volume and DetectorVolume (for gen1/3 and gen2 use cases)
std::pair<Volume, std::shared_ptr<Experimental::DetectorVolume>>
convertVolumeDetectorVolume(
    const GeometryContext& context, const GeoShape& shape,
    const std::string& name, const GeoTrf::Transform3D& transform,
    const std::vector<GeoModelSensitiveSurface>& sensitives);

}  // namespace Acts::GeoModel
