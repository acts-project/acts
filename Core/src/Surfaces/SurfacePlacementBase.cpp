// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfacePlacementBase.hpp"

#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

const Surface* SurfacePlacementBase::surface() const {
  return m_surface;
}

void SurfacePlacementBase::assignSurface(
    const std::shared_ptr<Surface>& surface) {
  m_surface = surface.get();

  const auto* placement = surface->surfacePlacement();

  // If the surface has no placement, assign this one to it
  if (placement == nullptr) {
    surface->assignSurfacePlacement(shared_from_this());
    return;
  }

  // If we get here, the surface has a placement, and it is not this one,
  // it's an error
  if (placement != this) {
    throw std::runtime_error(
        "SurfacePlacementBase: surface already has a placement. The placement "
        "cannot be reassigned");
  }
}

}  // namespace Acts