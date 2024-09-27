// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <tuple>
#include <utility>

#include "detray/core/detector.hpp"

namespace Acts::TracccPlugin {

/// @brief Creates a map from Acts geometry ID (value) to detray barcode.
/// @param detector the detray detector.
/// @return A map (key = geometry ID value, value = detray geometry barcode).
template <typename metadata_t, typename container_t>
inline std::map<Acts::GeometryIdentifier, detray::geometry::barcode>
createBarcodeMap(const detray::detector<metadata_t, container_t>& detector) {
  // Construct a map from Acts surface identifiers to Detray barcodes.
  std::map<Acts::GeometryIdentifier, detray::geometry::barcode> barcodeMap;
  for (const auto& surface : detector.surfaces()) {
    barcodeMap[Acts::GeometryIdentifier(surface.source)] = surface.barcode();
  }
  return barcodeMap;
}

}  // namespace Acts::TracccPlugin
