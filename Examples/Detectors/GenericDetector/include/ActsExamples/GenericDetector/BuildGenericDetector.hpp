// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingVolumeBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/GenericDetector/LayerBuilder.hpp"
#include "ActsExamples/GenericDetector/ProtoLayerCreator.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
class TrackingGeometry;
class HomogeneousSurfaceMaterial;
class IMaterialDecorator;
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples::Generic {

/// Global method to build the generic tracking geometry
///
/// @tparam detector_element_t is the actual type of the detector
/// element, each derivative of a GenericDetectorElement can be used
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorElementFactory is the factory for the detector elements
/// @param matDecorator is an optional decorator for the material
/// @param level is the detector building level
///          0 - pixel barrel only
///          1 - pixel detector only
///          2 - full barrel only
///          3 - full detector (without stereo modules)
/// @param matDecorator is the source for material decoration
/// @param protoMaterial is a flag to steer proto material to be loaded
/// @param surfaceLLevel is the surface building logging level
/// @param layerLLevel is the layer building logging level
/// @param volumeLLevel is the volume building logging level
/// return a unique vector to the tracking geometry
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const Acts::GeometryContext& gctxIn,
    const ProtoLayerCreator::DetectorElementFactory& detectorElementFactory,
    std::size_t level,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator = nullptr,
    bool protoMaterial = false,
    Acts::Logging::Level surfaceLLevel = Acts::Logging::INFO,
    Acts::Logging::Level layerLLevel = Acts::Logging::INFO,
    Acts::Logging::Level volumeLLevel = Acts::Logging::INFO);

}  // namespace ActsExamples::Generic
