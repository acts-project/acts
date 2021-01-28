// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include <array>
#include <memory>
#include <vector>

namespace ActsExamples {
namespace Telescope {

/// The telescope detector surface type
enum TelescopeSurfaceType : unsigned int {
  Plane = 0,
  Disc = 1,
};

/// Global method to build the telescope tracking geometry
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param positions is the offset w of different layers in the longitudinal
/// direction
/// @param offsets is the offset (u, v) of the layers in the transverse plane
/// @param bounds is the surface bound parameters, i.e. halfX and halfY if plane
/// surface, and minR and maxR if disc surface
/// @param thickness is the material thickness of each layer
/// @param surfaceType is the detector surface type
/// @param binValue indicates which axis the detector surface normals are
/// parallel to
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename TelescopeDetectorElement::ContextType& gctx,
    std::vector<std::shared_ptr<TelescopeDetectorElement>>& detectorStore,
    const std::vector<double>& positions, const std::array<double, 2>& offsets,
    const std::array<double, 2>& bounds, double thickness,
    const TelescopeSurfaceType& surfaceType,
    const Acts::BinningValue& binValue = Acts::BinningValue::binZ);

}  // end of namespace Telescope
}  // end of namespace ActsExamples
