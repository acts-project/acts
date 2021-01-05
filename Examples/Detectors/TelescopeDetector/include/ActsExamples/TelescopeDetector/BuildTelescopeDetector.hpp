// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
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

/// Global method to build the telescope tracking geometry
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param positions is the offset w of different layers in the longitudinal
/// direction
/// @param offsets is the offset (u, v) of the layers in the transverse plane
/// @param pSize is the plane size
/// @param thickness is the material thickness of each layer
/// @param binValue indicates which axis the planes normals are parallel to
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename TelescopeDetectorElement::ContextType& gctx,
    std::vector<std::shared_ptr<TelescopeDetectorElement>>& detectorStore,
    const std::vector<double>& positions, std::array<double, 2> offsets,
    std::array<double, 2> pSize, double thickness,
    Acts::BinningValue binValue = Acts::BinningValue::binZ);

}  // end of namespace Telescope
}  // end of namespace ActsExamples
