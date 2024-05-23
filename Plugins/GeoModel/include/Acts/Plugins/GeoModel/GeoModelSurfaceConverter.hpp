// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"

class GeoFullPhysVol;

namespace Acts {
/// Collect the sensitive surface & detector element
using GeoModelSensitiveSurface =
    std::tuple<std::shared_ptr<GeoModelDetectorElement>,
               std::shared_ptr<Surface>>;

namespace GeoModelSurfaceConverter {

// Options struct for the conversion
struct Options {
  Surface::SurfaceType cylinderTargetType = Surface::SurfaceType::Cylinder;
};

/// @brief conversion to sensitive surface
///
/// @param geoPhysVol the geoPhysVol to convert
/// @param options the conversion options
///
/// @note the conversion will identify the "thin" dimension of the shape
/// to allow a representation as a surface, if necessary the transform is
/// adjusted, such that the surface bounds are correctly defined.
/// This adjustment is attempted without flipping of axes as much as possible,
/// i.e. cyclic permutation of the axes is preferred over flipping.
///
///
/// @return a detector element and a surface
GeoModelSensitiveSurface convertToSensitiveSurface(
    const GeoFullPhysVol& geoPhysVol, const Options& options = {});
}  // namespace GeoModelSurfaceConverter
}  // namespace Acts
