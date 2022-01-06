// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>
#include <string>
#include <vector>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class GeometricExtent;
class DetectorVolume;
class InternalBlueprint;
class CylinderVolumeBounds;

struct CylindricalVolumeHelper {

/// Static helper method to build cylindrical bounds from a geometric extent
///
/// @param extent the geometric extent to build these bounds
/// @return a cylindrical volume bound object
static std::unique_ptr<CylinderVolumeBounds> buildBounds(
    const Acts::GeometricExtent& extent) noexcept(false);

/// Static helper method to build a transform from a geometric exent
///
/// @param extent the geometric extent to build the transform from
/// @return a transform to position the volume is global space
static Transform3 buildTransform(const Acts::GeometricExtent& extent) noexcept(
    false);

/// Static helper method to build a single volume
///
/// @param gctx the geometric context under which this is built
/// @param iBlueprint the internal blue print that is to be used for building
/// @param vExtent an external volume extent
/// @param name the name of the volume
///
/// @note this will throw an exception if the internal structure
///       is not confined by the external constraint (vExtent)
///
/// @return a newly built (shared) detector volume
static std::shared_ptr<DetectorVolume> buildVolume(
    const GeometryContext& gctx, const InternalBlueprint& iBlueprint,
    const GeometricExtent& vExtent, const std::string& name) noexcept(false);

/// Static helper method to build an ordered, sequenced list of volumes
///
/// @param gctx the geometric context under which this is built
/// @param iBlueprints the internal blue prints for the structured volumes
/// @param vExtent an external volume extent
/// @param keepValues the harmonized binning values, identical to all volumes
/// @param adaptValue the changing value for the volume creation
/// @param name the name of the volume
/// @param volSuffix a string identifying structured sub volumes
/// @param gapSurffix a string indentifying gap volumes
/// @param logLevel a buildint log level for screen output
///
/// @note this will throw an exception if any of the internal structurse
///       are not confined by the external constraint (vExtent)
/// @note this will also throw an exception if no blueprints are provided
///
/// @return a newly built (shared) detector volume
static std::vector<std::shared_ptr<DetectorVolume>> buildVolumes(
    const GeometryContext& gctx,
    const std::vector<InternalBlueprint>& iBlueprints,
    const GeometricExtent& vExtent, const std::vector<BinningValue>& keepValues,
    const Acts::BinningValue& adaptValue, const std::string& name,
    const std::string& volSuffix, const std::string& gapSuffix,
    Acts::Logging::Level logLevel) noexcept(false);

/// Static helper method to check the volume and blueprint extents
///
/// @param external is the external volume constraint
/// @param iBlueprints are the internal structure blueprints
///
/// @note this will throw and exception if any of the internal structures
///       are not contained by the external one
static void checkExtents(
    const GeometricExtent& external,
    const std::vector<InternalBlueprint>& iBlueprints) noexcept(false);

// Static helper method to harmonize the internal layer extends
///
/// @param external is the external volume constraint
/// @param iBlueprints are the internal structure blueprints
/// @param binValues are the binning values which will be harmonized/set
///
/// @note this will method will modify the extents of the internal blue prints
///       if necessary
///
/// @return a new extend for the container
static GeometricExtent harmonizeExtents(
    const GeometricExtent& external,
    std::vector<InternalBlueprint>& iBlueprints,
    const std::vector<BinningValue>& binValues = s_binningValues);

};

}  // namespace Acts
