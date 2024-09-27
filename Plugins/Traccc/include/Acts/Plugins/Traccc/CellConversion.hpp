// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Traccc/BarcodeMap.hpp"
#include "Acts/Plugins/Traccc/DigitizationConfig.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <tuple>
#include <utility>

#include "detray/core/detector.hpp"
#include "traccc/edm/cell.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"

namespace Acts::TracccPlugin {

/// @brief Converts a "geometry ID -> traccc cells" map to traccc cells and modules.
/// @param mr The memory resource to use.
/// @param cellsMap A map from Acts geometry ID value to traccc cells.
/// @param geom The traccc geometry.
/// @param dconfig The traccc digitization configuration.
/// @param barcode_map A map from Acts geometry ID value to detray barcode.
/// @return A tuple containing the traccc cells (first item) and traccc modules (second item).
std::tuple<traccc::cell_collection_types::host,
           traccc::silicon_detector_description::host>
createCellsAndModules(
    vecmem::memory_resource& mr,
    const std::map<Acts::GeometryIdentifier, std::vector<traccc::cell>>&
        cellsMap,
    const std::map<Acts::GeometryIdentifier, traccc::transform3>& geom,
    const DigitizationConfig& dconfig,
    const std::map<Acts::GeometryIdentifier, detray::geometry::barcode>&
        barcodeMap,
    std::unique_ptr<const Acts::Logger> logger);

}  // namespace Acts::TracccPlugin
