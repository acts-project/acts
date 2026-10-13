// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <vector>

namespace ActsExamples {

/// One module's entry in the layer mapping file.
struct GbtsLayerInfo {
  /// GBTS layer id
  Acts::Experimental::GbtsExperimentLayerId layerId{};
  /// whether the layer is a barrel or an endcap layer
  Acts::Experimental::GbtsLayerType type{};
  /// sensor technology of the layer
  Acts::Experimental::GbtsLayerTechnology technology{};
};

/// conversion between ACTS labelling of volume, layer and modules to that used
/// by GBTS: the entry of a surface is the one of its module or, without one,
/// the one of its whole layer
using GbtsLayerMap = Acts::GeometryHierarchyMap<GbtsLayerInfo>;

/// make the map between ACTS geometry ID's and GBTS geometry ID's
/// @param layerMappingFile The layer mapping file to read
/// @return The entry of every surface of the file
GbtsLayerMap makeGbtsLayerMap(const std::filesystem::path &layerMappingFile);

/// makes the geometry objects used by GBTS that correspond to the objects in
/// the connection table for ease these are sometimes called "logical layers"
/// @param trackingGeometry The detector the layers are made of
/// @param layerMap The GBTS layer of the surfaces
/// @param gctx The geometry context
/// @param useLayerTilt Whether to give inclined layers the slope of their
///                     modules, along which the eta bins are cut
/// @param fillModuleCsv Whether to write every module to ACTS_modules.csv
/// @param logger The logger the surfaces on none of the layers go to
/// @return One description per GBTS layer
std::vector<Acts::Experimental::GbtsLayerDescription> makeGbtsLayerDescriptions(
    const Acts::TrackingGeometry &trackingGeometry,
    const GbtsLayerMap &layerMap, const Acts::GeometryContext &gctx,
    bool useLayerTilt, bool fillModuleCsv, const Acts::Logger &logger);

}  // namespace ActsExamples
