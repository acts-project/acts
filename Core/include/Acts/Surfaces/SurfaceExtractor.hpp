// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"

#include <vector>

namespace Acts {

/// @brief selector for extracting surfaces, this is done without
/// having a requirement on geometry indentifier
struct SurfaceExtractor {
  bool material = true;
  bool sensitive = false;
  bool all = false;

  std::vector<const Surface*> extractedSurfaces;

  /// @brief selector for finding surface
  ///
  /// @param surface is the surface to be checked/visited
  void operator()(const Surface* surface) {
    if (all || (material && surface->surfaceMaterial() != nullptr)) {
      extractedSurfaces.push_back(surface);
      return;
    }
    if (sensitive && surface->associatedDetectorElement() != nullptr) {
      extractedSurfaces.push_back(surface);
      return;
    }
  }
};

}  // namespace Acts
