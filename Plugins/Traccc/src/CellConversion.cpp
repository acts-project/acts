// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts include(s)
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Plugin include(s)
#include "Acts/Plugins/Traccc/BarcodeMap.hpp"
#include "Acts/Plugins/Traccc/DigitizationConfig.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"

// Traccc include(s)
#include "traccc/edm/cell.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"

// VecMem include(s)
#include "vecmem/memory/memory_resource.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <tuple>
#include <utility>

namespace {

/// Comparator used for sorting cells. This sorting is one of the assumptions
/// made in the clusterization algorithm
struct CellOrder {
  bool operator()(const traccc::cell& lhs, const traccc::cell& rhs) const {
    if (lhs.module_link != rhs.module_link) {
      return lhs.module_link < rhs.module_link;
    } else if (lhs.channel1 != rhs.channel1) {
      return (lhs.channel1 < rhs.channel1);
    } else {
      return (lhs.channel0 < rhs.channel0);
    }
  }
};  // struct CellOrder

}  // namespace

namespace Acts::TracccPlugin {

std::tuple<traccc::cell_collection_types::host,
           traccc::silicon_detector_description::host>
createCellsAndModules(
    vecmem::memory_resource& mr,
    std::map<Acts::GeometryIdentifier::Value, std::vector<traccc::cell>>&
        cellsMap,
    std::optional<std::reference_wrapper<const traccc::geometry>> geom,
    const DigitizationConfig& dconfig,
    const std::map<Acts::GeometryIdentifier::Value, detray::geometry::barcode>&
        barcodeMap) {
  traccc::cell_collection_types::host cells(&mr);
  traccc::silicon_detector_description::host modules(mr);

  // Sort the cells.
  for (auto& [_, mapCells] : cellsMap) {
    std::ranges::sort(mapCells, CellOrder());
  }

  modules.resize(cellsMap.size());

  std::size_t module_idx = 0;

  // Fill the output containers with the ordered cells and modules.
  for (const auto& [originalGeometryID, mapCells] : cellsMap) {
    // Modify the geometry ID of the module if a barcode map is
    // provided.
    Acts::GeometryIdentifier::Value geometryID =
        barcodeMap.at(originalGeometryID).value();

    // Add the module and its cells to the output.
    modules.geometry_id().at(module_idx) =
        detray::geometry::barcode{geometryID};
    modules.acts_geometry_id().at(module_idx) = originalGeometryID;
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
    if (geom) {
      const traccc::geometry& g = (*geom).get();

      // Check if the module ID is known.
      if (g.contains(originalGeometryID)) {
        throw std::runtime_error("Could not find placement for geometry ID " +
                                 std::to_string(originalGeometryID));
      }

      // Set the value on the module description.
      modules.placement().at(module_idx) = g.at(originalGeometryID);
    }

    const auto binning_data = digi_it->segmentation.binningData();
    modules.reference_x().at(module_idx) = binning_data.at(0).min;
    modules.reference_y().at(module_idx) = binning_data.at(1).min;
    modules.pitch_x().at(module_idx) = binning_data.at(0).step;
    modules.pitch_y().at(module_idx) = binning_data.at(1).step;
    modules.dimensions().back() = digi_it->dimensions;

    module_idx++;

    for (auto& cell : mapCells) {
      cells.push_back(cell);
      // Set the module link.
      cells.back().module_link = modules.size() - 1;
    }
  }
  return std::make_tuple(std::move(cells), std::move(modules));
}

}  // namespace Acts::TracccPlugin
