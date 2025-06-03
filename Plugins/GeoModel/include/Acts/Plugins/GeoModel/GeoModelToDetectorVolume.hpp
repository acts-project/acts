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
#include "Acts/Utilities/BoundFactory.hpp"

#include "GeoModelKernel/GeoDefinitions.h"
class GeoShape;

namespace Acts::GeoModel {

/// @brief Converts a GeoShape into a bounded volume. For the supported shape types and the
///        chosen strategie please consult the cpp file. May throw an exception
/// f the GeoShape is not yet considered.
/// @param trf: Transform to align position in the volume in space
/// @param shape: Pointer to the GeoShape from which the VolumeBounds are translated
/// @param boundFactory: Reference to the bound factory to avoid multiple instances of
///                      equivalent bound parameters
/// @return A shared pointer initialized with the new volume
std::shared_ptr<Volume> convertVolume(const Transform3& trf,
                                      const GeoShape* shape,
                                      VolumeBoundFactory& boundFactory);

/// @brief Converts a simple Volume into a Gen-2 DetectorVolume with associated sensitive surfaces inside
/// @param context: GeometryContext to align the volume needed during the construction phase of the volume
/// @param vol: Non-const reference to the volume to convert
/// @param name: Name of the constructed Gen-2 volume
/// @param sensitives: List of sensitive surfaces to be put inside the detector volume.
/// @return A shared pointer initialized with the new Gen-2 volume
std::shared_ptr<Experimental::DetectorVolume> convertDetectorVolume(
    const GeometryContext& context, Volume& vol, const std::string& name,
    const std::vector<std::shared_ptr<Surface>>& sensitives);

}  // namespace Acts::GeoModel
