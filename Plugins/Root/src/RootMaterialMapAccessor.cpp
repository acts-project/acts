// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMaterialMapAccessor.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IGridSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

void Acts::RootMaterialMapAccessor::write(TFile& rFile,
                                          const GeometryContext& gctx,
                                          const Surface& surface) const {
  auto geoID = surface.geometryId();
  auto surfaceMaterial = surface.surfaceMaterial();
  if (surfaceMaterial == nullptr) {
    return;
  }

  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(surfaceMaterial);
  if (homogeneousMaterial != nullptr) {
    
    
  }

  auto binnedMaterial =
      dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
  if (binnedMaterial != nullptr) {
    // writeBinnedMaterial(rFile, gctx, geoID, *binnedMaterial);
  }

  return;
}
