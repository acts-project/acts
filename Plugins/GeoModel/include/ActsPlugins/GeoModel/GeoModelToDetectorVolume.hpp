// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Utilities/BoundFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

#include "GeoModelKernel/GeoDefinitions.h"
class GeoShape;

namespace ActsPlugins::GeoModel {

/// @addtogroup geomodel_plugin
/// @{

/// @brief Calculates the absolute volume position w.r.t. the world node
///        This is only possible, if the volume is not shared in multiple
///        branches of the GeoModelTree.
/// @param physVol: Reference to the physical volume from which the position
///                 shall be calculated
Acts::Transform3 volumePosInSpace(const PVConstLink& physVol);
/// @brief Converts a GeoShape into a bounded volume. For the supported shape types and the
///        chosen strategie please consult the cpp file. May throw an exception
/// f the GeoShape is not yet considered.
/// @param trf: Transform to align position in the volume in space
/// @param shape: Pointer to the GeoShape from which the VolumeBounds are translated
/// @param boundFactory: Reference to the bound factory to avoid multiple instances of
///                      equivalent bound parameters
/// @return A shared pointer initialized with the new volume
std::shared_ptr<Acts::Volume> convertVolume(
    const Acts::Transform3& trf, const GeoShape* shape,
    Acts::VolumeBoundFactory& boundFactory);

/// @}

}  // namespace ActsPlugins::GeoModel
