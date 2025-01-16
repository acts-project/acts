// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <map>

#include "detray/core/detector.hpp"
#include "traccc/definitions/primitives.hpp"

namespace Acts::TracccPlugin {

/// @brief Creates a map from Acts geometry ID (value) to detray barcode.
/// @param detector the detray detector.
/// @return A map (key = geometry ID value, value = detray geometry barcode).
template <typename metadata_t, typename container_t>
std::map<Acts::GeometryIdentifier, traccc::transform3> createSurfaceMap(
    const detray::detector<metadata_t, container_t>& detector) {
  std::map<Acts::GeometryIdentifier, traccc::transform3> maps;

  const typename detray::detector<metadata_t, container_t>::geometry_context
      ctx0{};

  for (const auto& sf_desc : detector.surfaces()) {
    const detray::tracking_surface sf{detector, sf_desc.barcode()};

    maps.insert({Acts::GeometryIdentifier(sf_desc.source), sf.transform(ctx0)});
  }

  return maps;
}

}  // namespace Acts::TracccPlugin
