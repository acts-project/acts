// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/Root/RootMaterialMapAccessor.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <memory>
#include <tuple>
#include <vector>

#include "TFile.h"

using namespace Acts;

using IdentifiedMaterial =
    std::tuple<GeometryIdentifier, std::shared_ptr<ISurfaceMaterial>>;

std::vector<IdentifiedMaterial> createHomogeneousSurfaceMaterial() {
  std::size_t nMaterials = 100;

  std::vector<IdentifiedMaterial> homogeneousMaterials;
  homogeneousMaterials.reserve(nMaterials);
  for (std::size_t i = 0; i < nMaterials; ++i) {
    // construct the material properties from arguments
    Material mat = Material::fromMolarDensity(
        1. + i * 0.5, 2. + i * 0.5, 3. + i * 0.5, 4. + i * 0.5, 5. + i * 0.5);
    MaterialSlab mp(mat, 0.1);
    auto hMaterial = std::make_shared<HomogeneousSurfaceMaterial>(mp);
    auto geoID = GeometryIdentifier().withVolume(1).withSensitive(i + 1);
    homogeneousMaterials.push_back({geoID, hMaterial});
  }
  return homogeneousMaterials;
}

std::vector<IdentifiedMaterial> createBinnedSurfaceMaterial() {
  std::size_t nMaterials = 100;

  std::vector<IdentifiedMaterial> binnedMaterials;
  binnedMaterials.reserve(nMaterials);
  for (std::size_t i = 0; i < nMaterials; ++i) {
    // construct the material properties from arguments

    BinUtility xyBinning(100, -1., 1., open, AxisDirection::AxisX);
    xyBinning += BinUtility(50, -3., 3., open, AxisDirection::AxisY);

    std::vector<std::vector<MaterialSlab>> materialMatrix;
    for (std::size_t j = 0; j < xyBinning.bins(1); ++j) {
      std::vector<MaterialSlab> materialRow;
      for (std::size_t k = 0; k < xyBinning.bins(0); ++k) {
        // Create a material slab with some arbitrary properties
        Material mat = Material::fromMolarDensity(
            i + j * 1. + k * 0.5, i + j * 2 + k * 0.5, i + j * 3. + k * 0.5,
            i + j * 4. + k * 0.5, i + j * 5. + k * 0.5);
        MaterialSlab mp(mat, 0.1);
        materialRow.push_back(mp);
      }
      materialMatrix.push_back(materialRow);
    }
    auto binnedMaterial =
        std::make_shared<BinnedSurfaceMaterial>(xyBinning, materialMatrix);
    auto geoID = GeometryIdentifier().withVolume(2).withSensitive(i + 1);
    binnedMaterials.push_back({geoID, binnedMaterial});
  }
  return binnedMaterials;
}

BOOST_AUTO_TEST_SUITE(RootMaterialMapAccessorTests)

BOOST_AUTO_TEST_CASE(RootMaterialMapAccessorHomogeneousReadWrite) {
  auto surfaceMaterials = createHomogeneousSurfaceMaterial();

  auto rFile =
      TFile::Open("RootMaterialMapAccessorHomogeneousTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  // Create the accessor
  Acts::RootMaterialMapAccessor::Config cfg;
  Acts::RootMaterialMapAccessor accessor(cfg);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessor.write(*rFile, geoID, *sMaterial);
  }

  rFile->Write();
  rFile->Close();

  // Let's read it back
  auto iFile =
      TFile::Open("RootMaterialMapAccessorHomogeneousTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);

  auto [surfaceMapsRead, volumeMapsRead] = accessor.read(*iFile);
  BOOST_REQUIRE_EQUAL(surfaceMapsRead.size(), surfaceMaterials.size());
  BOOST_REQUIRE_EQUAL(volumeMapsRead.size(), 0);

  Vector3 accessorPosition(0., 0., 0.);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    auto it = surfaceMapsRead.find(geoID);
    BOOST_REQUIRE(it != surfaceMapsRead.end());
    const auto& readMaterial = it->second;
    BOOST_REQUIRE(readMaterial != nullptr);
    const auto* hMaterial =
        dynamic_cast<const HomogeneousSurfaceMaterial*>(readMaterial.get());
    BOOST_REQUIRE(hMaterial != nullptr);
    BOOST_CHECK_CLOSE(hMaterial->materialSlab(accessorPosition).material().X0(),
                      sMaterial->materialSlab(accessorPosition).material().X0(), 1e-6);
    BOOST_CHECK_CLOSE(hMaterial->materialSlab(accessorPosition).material().L0(),
                      sMaterial->materialSlab(accessorPosition).material().L0(), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(RootMaterialMapAccessorBinnrfReadWrite) {
  auto surfaceMaterials = createBinnedSurfaceMaterial();

  auto rFile =
      TFile::Open("RootMaterialMapAccessorBinnedTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  // Create the accessor
  Acts::RootMaterialMapAccessor::Config cfg;
  Acts::RootMaterialMapAccessor accessor(cfg);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessor.write(*rFile, geoID, *sMaterial);
  }

  rFile->Write();
  rFile->Close();

  // Create the accessor - compressed writing
  Acts::RootMaterialMapAccessor::Config cfgIndexed{true};
  Acts::RootMaterialMapAccessor accessorIndexed(cfgIndexed);

  rFile =
      TFile::Open("RootMaterialMapAccessorBinnedIndexedTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessorIndexed.write(*rFile, geoID, *sMaterial);
  }

  rFile->Write();
  rFile->Close();
}

BOOST_AUTO_TEST_SUITE_END()
