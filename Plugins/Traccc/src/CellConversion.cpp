// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
#include "vecmem/memory/memory_resource.hpp"

namespace Acts::TracccPlugin {

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
    std::unique_ptr<const Acts::Logger> _logger) {
  ACTS_LOCAL_LOGGER(std::move(_logger));

  traccc::cell_collection_types::host cells(&mr);
  traccc::silicon_detector_description::host modules(mr);

  modules.resize(cellsMap.size());

  std::size_t module_idx = 0;

  // Fill the output containers with the ordered cells and modules.
  for (const auto& [originalGeometryID, mapCells] : cellsMap) {
    // Modify the geometry ID of the module if a barcode map is
    // provided.
    auto it = barcodeMap.find(originalGeometryID);

    if (it == barcodeMap.cend()) {
      ACTS_WARNING("Geometry element " << originalGeometryID
                                       << " was not found in the barcode map.");
      continue;
    }

    Acts::GeometryIdentifier::Value geometryID = it->second.value();

    // Add the module and its cells to the output.
    modules.geometry_id().at(module_idx) =
        detray::geometry::barcode{geometryID};
    modules.acts_geometry_id().at(module_idx) = originalGeometryID.value();
    modules.threshold().at(module_idx) = 0.f;

    const DigitizationConfig::Iterator digi_it =
        dconfig.find(originalGeometryID);

    if (digi_it == dconfig.end()) {
      std::ostringstream msg;
      msg << "Could not find digitization config for geometry ID: "
          << originalGeometryID;
      throw std::runtime_error(msg.str());
    }

    // Find/set the 3D position of the detector module.
    if (!geom.contains(originalGeometryID.value())) {
      throw std::runtime_error("Could not find placement for geometry ID " +
                               std::to_string(originalGeometryID.value()));
    }

    // Set the value on the module description.
    modules.placement().at(module_idx) = geom.at(originalGeometryID.value());

    const auto binning_data = digi_it->segmentation.binningData();

    if (!binning_data.empty()) {
      modules.reference_x().at(module_idx) = binning_data.at(0).min;
      modules.pitch_x().at(module_idx) = binning_data.at(0).step;
    } else {
      modules.reference_x().at(module_idx) = 0.f;
      modules.pitch_x().at(module_idx) = 1.f;
    }

    if (binning_data.size() >= 2) {
      modules.reference_y().at(module_idx) = binning_data.at(1).min;
      modules.pitch_y().at(module_idx) = binning_data.at(1).step;
    } else {
      modules.reference_y().at(module_idx) = 0.f;
      modules.pitch_y().at(module_idx) = 1.f;
    }

    modules.dimensions().at(module_idx) = digi_it->dimensions;

    for (auto& cell : mapCells) {
      cells.push_back(cell);
      // std::cout << cell.channel0 << ", " << cell.channel1 << " @ " << module_idx << std::endl;
      // Set the module link.
      cells.back().module_link = module_idx;
    }

    module_idx++;
  }
  return std::make_tuple(std::move(cells), std::move(modules));
}

}  // namespace Acts::TracccPlugin
