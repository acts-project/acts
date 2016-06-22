// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_BUILDGENERICDETECTOR_H
#define ACTS_BUILDGENERICDETECTOR_H 1

#include <vector>
#include <memory>
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {
class TrackingGeometry;

/// Global method to build the generic tracking geometry
/// @param lvl is the debug logging level
/// @param version is the detector version
std::unique_ptr<const Acts::TrackingGeometry> 
trackingGeometry(Acts::Logging::Level lvl = Acts::Logging::INFO, size_t version=0);

/// Helper method for positioning
/// @param radius is the cylinder radius
/// @param tilt is the nominal layer radius
/// @param zStagger is the radial staggering along z
/// @param lModule is the module length (longitudinal)
/// @param lOverlap is the overlap of the modules (longitudinal)
/// @binningSchema is the way the bins are laid out rphi x z
std::vector<  Acts::Vector3D >
modulePositionsCylinder(double radius,
                        double zStagger,
                        double lModule,
                        double lOverlap,
                        const std::pair<int,int>& binningSchema);

}

#endif  // ACTS_BUILDGENERICDETECTOR_H
