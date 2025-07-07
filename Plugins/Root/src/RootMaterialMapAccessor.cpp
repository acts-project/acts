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
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

void Acts::RootMaterialMapAccessor::write(TFile& rFile,
                                          const GeometryContext& gctx,
                                          const Surface& surface)  {
  auto geoID = surface.geometryId();
  auto surfaceMaterial = surface.surfaceMaterial();
  if (surfaceMaterial == nullptr) {
    return;
  }


  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(surfaceMaterial);
  if (homogeneousMaterial != nullptr) {
    writeHomogeneousMaterial(*homogeneousMaterial);
    m_hGeoIdPtr->push_back(geoID.value());
  }

  auto binnedMaterial =
      dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
  if (binnedMaterial != nullptr) {
    // writeBinnedMaterial(rFile, gctx, geoID, *binnedMaterial);
  }

  return;
}

void Acts::RootMaterialMapAccessor::connectForWrite(
    TTree& rTree, MaterialTreePayload& treePayload) {
  rTree.Branch("hGeoId", &treePayload.hGeoIdPtr);
  rTree.Branch("ht", &treePayload.htPtr);
  rTree.Branch("hX0", &treePayload.hX0Ptr);
  rTree.Branch("hL0", &treePayload.hL0Ptr);
  rTree.Branch("hA", &treePayload.hAPtr);
  rTree.Branch("hZ", &treePayload.hZPtr);
  rTree.Branch("hRho", &treePayload.hRhoPtr);
}

void Acts::RootMaterialMapAccessor::connectForRead(
    const TTree& rTree, MaterialTreePayload& treePayload) {
  rTree.SetBranchAddress("hGeoId", &treePayload.hGeoIdPtr);
  rTree.SetBranchAddress("ht", &treePayload.htPtr);
  rTree.SetBranchAddress("hX0", &treePayload.hX0Ptr);
  rTree.SetBranchAddress("hL0", &treePayload.hL0Ptr);
  rTree.SetBranchAddress("hA", &treePayload.hAPtr);
  rTree.SetBranchAddress("hZ", &treePayload.hZPtr);
  rTree.SetBranchAddress("hRho", &treePayload.hRhoPtr);
}

void Acts::RootMaterialMapAccessor::writeHomogeneousMaterial(
    const HomogeneousSurfaceMaterial& homogeneousMaterial) {
  if (m_hTree == nullptr) {
    m_hTree = std::make_unique<TTree>("HomogeneousMaterial",
                                      "Homogeneous Material Tree");


  }


}
