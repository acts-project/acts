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
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace ActsPlugins;

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

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(RootMaterialMapIoHomogeneousReadWrite) {
  auto surfaceMaterials = createHomogeneousSurfaceMaterial();

  auto rFile =
      TFile::Open("RootMaterialMapIoHomogeneousTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  // Create the accessor
  RootMaterialMapIo::Config cfg;
  RootMaterialMapIo accessor(cfg);
  RootMaterialMapIo::Options options;

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessor.write(*rFile, geoID, *sMaterial, options);
  }

  rFile->Write();
  rFile->Close();

  // Let's read it back
  auto iFile = TFile::Open("RootMaterialMapIoHomogeneousTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);

  auto [surfaceMapsRead, volumeMapsRead] = accessor.read(*iFile, options);
  BOOST_REQUIRE_EQUAL(surfaceMapsRead.size(), surfaceMaterials.size());
  BOOST_REQUIRE_EQUAL(volumeMapsRead.size(), 0);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    auto it = surfaceMapsRead.find(geoID);
    BOOST_REQUIRE(it != surfaceMapsRead.end());
    const auto& readMaterial = it->second;
    BOOST_REQUIRE(readMaterial != nullptr);
    const auto* hMaterial =
        dynamic_cast<const HomogeneousSurfaceMaterial*>(readMaterial.get());
    BOOST_REQUIRE(hMaterial != nullptr);
    BOOST_CHECK_CLOSE(hMaterial->materialSlab().material().X0(),
                      sMaterial->materialSlab(Vector2{0., 0.}).material().X0(),
                      1e-6);
    BOOST_CHECK_CLOSE(hMaterial->materialSlab().material().L0(),
                      sMaterial->materialSlab(Vector2{0., 0.}).material().L0(),
                      1e-6);
  }
}

// Legacy BinnedSurfaceMaterial files (both the raw-histogram and the
// indexed on-disk shape) must be readable, but are always upgraded to a
// GridSurfaceMaterial on read - RootMaterialMapIo no longer reconstructs
// BinnedSurfaceMaterial at all.
BOOST_AUTO_TEST_CASE(RootMaterialMapIoBinnedReadAsGrid) {
  auto surfaceMaterials = createBinnedSurfaceMaterial();

  auto rFile = TFile::Open("RootMaterialMapIoBinnedTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  // Create the accessor
  RootMaterialMapIo::Config cfg;
  RootMaterialMapIo accessor(cfg);
  RootMaterialMapIo::Options options;

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessor.write(*rFile, geoID, *sMaterial, options);
  }

  rFile->Write();
  rFile->Close();

  // Let's read it back
  auto iFile = TFile::Open("RootMaterialMapIoBinnedTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);
  auto [surfaceMapsRead, volumeMapsRead] = accessor.read(*iFile, options);
  BOOST_REQUIRE_EQUAL(surfaceMapsRead.size(), surfaceMaterials.size());
  BOOST_REQUIRE_EQUAL(volumeMapsRead.size(), 0);

  // Bin centers of the 100 x [-1,1] / 50 x [-3,3] binning used above
  auto xCenter = [](std::size_t k) { return -1. + (k + 0.5) * 2. / 100.; };
  auto yCenter = [](std::size_t j) { return -3. + (j + 0.5) * 6. / 50.; };

  // Compare
  for (const auto& [refGeoID, refSMaterial] : surfaceMaterials) {
    auto binnedReferenceMaterial =
        dynamic_cast<const BinnedSurfaceMaterial*>(refSMaterial.get());
    BOOST_REQUIRE(binnedReferenceMaterial != nullptr);

    auto it = surfaceMapsRead.find(refGeoID);
    BOOST_REQUIRE(it != surfaceMapsRead.end());
    const auto& readMaterial = it->second;
    BOOST_REQUIRE(readMaterial != nullptr);
    const auto* gridMaterial =
        dynamic_cast<const GridSurfaceMaterial*>(readMaterial.get());
    BOOST_REQUIRE(gridMaterial != nullptr);

    const auto& referenceMaterialMatrix =
        binnedReferenceMaterial->fullMaterial();
    BOOST_REQUIRE_EQUAL(referenceMaterialMatrix.size(), 50u);
    for (std::size_t j = 0; j < 50; ++j) {
      BOOST_REQUIRE_EQUAL(referenceMaterialMatrix[j].size(), 100u);
      for (std::size_t k = 0; k < 100; ++k) {
        const auto& refMat = referenceMaterialMatrix[j][k];
        const auto& mat =
            gridMaterial->materialSlab(Vector2{xCenter(k), yCenter(j)});
        BOOST_CHECK_CLOSE(mat.material().X0(), refMat.material().X0(), 1e-4);
        BOOST_CHECK_CLOSE(mat.material().L0(), refMat.material().L0(), 1e-4);
        BOOST_CHECK_CLOSE(mat.material().Ar(), refMat.material().Ar(), 1e-4);
        BOOST_CHECK_CLOSE(mat.material().Z(), refMat.material().Z(), 1e-4);
        BOOST_CHECK_CLOSE(mat.thickness(), refMat.thickness(), 1e-4);
      }
    }
  }

  // Same, but written through the indexed BinnedSurfaceMaterial mode
  RootMaterialMapIo::Config cfgIndexed;
  RootMaterialMapIo accessorIndexed(cfgIndexed);

  RootMaterialMapIo::Options optionsIndexed;
  optionsIndexed.indexedMaterial = true;

  rFile = TFile::Open("RootMaterialMapIoBinnedIndexedTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    accessorIndexed.write(*rFile, geoID, *sMaterial, optionsIndexed);
  }

  rFile->Write();
  rFile->Close();

  // Let's read it back
  iFile = TFile::Open("RootMaterialMapIoBinnedIndexedTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);
  auto [surfaceMapsIndexedRead, volumeMapsIndexedRead] =
      accessorIndexed.read(*iFile, optionsIndexed);
  BOOST_REQUIRE_EQUAL(surfaceMapsIndexedRead.size(), surfaceMaterials.size());
  BOOST_REQUIRE_EQUAL(volumeMapsIndexedRead.size(), 0);

  for (const auto& [refGeoID, refSMaterial] : surfaceMaterials) {
    auto binnedReferenceMaterial =
        dynamic_cast<const BinnedSurfaceMaterial*>(refSMaterial.get());
    BOOST_REQUIRE(binnedReferenceMaterial != nullptr);
    auto it = surfaceMapsIndexedRead.find(refGeoID);
    BOOST_REQUIRE(it != surfaceMapsIndexedRead.end());
    const auto& readMaterial = it->second;
    BOOST_REQUIRE(readMaterial != nullptr);
    const auto* gridMaterial =
        dynamic_cast<const GridSurfaceMaterial*>(readMaterial.get());
    BOOST_REQUIRE(gridMaterial != nullptr);

    const auto& referenceMaterialMatrix =
        binnedReferenceMaterial->fullMaterial();
    for (std::size_t j = 0; j < 50; ++j) {
      for (std::size_t k = 0; k < 100; ++k) {
        const auto& refMat = referenceMaterialMatrix[j][k];
        const auto& mat =
            gridMaterial->materialSlab(Vector2{xCenter(k), yCenter(j)});
        BOOST_CHECK_CLOSE(mat.material().X0(), refMat.material().X0(), 1e-4);
        BOOST_CHECK_CLOSE(mat.thickness(), refMat.thickness(), 1e-4);
      }
    }
  }
}

// A natively written GridSurfaceMaterial (direct storage in memory) must
// round-trip through the file, reconstructed with the default (Indexed)
// storage kind.
BOOST_AUTO_TEST_CASE(RootMaterialMapIoGridDirectReadWrite) {
  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 4);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, -2., 2., 3);

  std::vector<std::vector<MaterialSlab>> payload;
  for (std::size_t i0 = 0; i0 < 4; ++i0) {
    std::vector<MaterialSlab> row;
    for (std::size_t i1 = 0; i1 < 3; ++i1) {
      Material mat = Material::fromMolarDensity(1. + i0, 2. + i1, 3. + i0,
                                                4. + i1, 5. + i0);
      row.emplace_back(mat, 1. + 0.1 * static_cast<double>(i0 + i1));
    }
    payload.push_back(row);
  }
  auto gridMaterial =
      GridSurfaceMaterialFactory::createDirect(*axisX, *axisY, payload);
  auto geoID = GeometryIdentifier().withVolume(3).withSensitive(1);

  auto rFile = TFile::Open("RootMaterialMapIoGridDirectTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  RootMaterialMapIo::Config cfg;
  RootMaterialMapIo accessor(cfg);
  RootMaterialMapIo::Options options;
  accessor.write(*rFile, geoID, *gridMaterial, options);
  rFile->Write();
  rFile->Close();

  auto iFile = TFile::Open("RootMaterialMapIoGridDirectTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);
  auto [surfaceMapsRead, volumeMapsRead] = accessor.read(*iFile, options);
  BOOST_REQUIRE_EQUAL(surfaceMapsRead.size(), 1u);

  auto it = surfaceMapsRead.find(geoID);
  BOOST_REQUIRE(it != surfaceMapsRead.end());
  const auto* readGridMaterial =
      dynamic_cast<const GridSurfaceMaterial*>(it->second.get());
  BOOST_REQUIRE(readGridMaterial != nullptr);
  // The persisted MaterialKind is authoritative, so this reads back as
  // Direct regardless of options.gridStorageKind (that only applies to
  // legacy BinnedSurfaceMaterial data, which carries no such marker)
  BOOST_CHECK(std::holds_alternative<GridSurfaceMaterial::Direct>(
      readGridMaterial->storage()));

  for (std::size_t i0 = 0; i0 < 4; ++i0) {
    double x = -1. + (i0 + 0.5) * 2. / 4.;
    for (std::size_t i1 = 0; i1 < 3; ++i1) {
      double y = -2. + (i1 + 0.5) * 4. / 3.;
      const auto& refMat = payload[i0][i1];
      const auto& mat = readGridMaterial->materialSlab(Vector2{x, y});
      BOOST_CHECK_CLOSE(mat.material().X0(), refMat.material().X0(), 1e-4);
      BOOST_CHECK_CLOSE(mat.thickness(), refMat.thickness(), 1e-4);
    }
  }
}

// Two GridSurfaceMaterial instances sharing the same GloballyIndexed material
// vector must have that vector written to the shared tree only once.
BOOST_AUTO_TEST_CASE(RootMaterialMapIoGridGloballyIndexedDedup) {
  auto sharedMaterial = std::make_shared<std::vector<MaterialSlab>>();
  sharedMaterial->emplace_back(Material::Vacuum(), 0.0);
  sharedMaterial->emplace_back(Material::fromMolarDensity(1., 2., 3., 4., 5.),
                               1.0);
  sharedMaterial->emplace_back(
      Material::fromMolarDensity(11., 12., 13., 14., 15.), 2.0);

  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 2., 2);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 1., 1);

  std::vector<std::vector<std::size_t>> indexPayload0 = {{1u}, {2u}};
  std::vector<std::vector<std::size_t>> indexPayload1 = {{2u}, {1u}};

  auto gridMaterial0 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axisX, *axisY, sharedMaterial, indexPayload0);
  auto gridMaterial1 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axisX, *axisY, sharedMaterial, indexPayload1);

  auto geoID0 = GeometryIdentifier().withVolume(4).withSensitive(1);
  auto geoID1 = GeometryIdentifier().withVolume(4).withSensitive(2);

  auto rFile =
      TFile::Open("RootMaterialMapIoGridGloballyIndexedTests.root", "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  RootMaterialMapIo::Config cfg;
  RootMaterialMapIo accessor(cfg);
  RootMaterialMapIo::Options options;
  accessor.write(*rFile, geoID0, *gridMaterial0, options);
  accessor.write(*rFile, geoID1, *gridMaterial1, options);
  rFile->Write();

  // The shared material vector must only have been written once: 3 rows,
  // not 6
  auto* tree =
      dynamic_cast<TTree*>(rFile->Get(options.indexedMaterialTreeName.c_str()));
  BOOST_REQUIRE(tree != nullptr);
  BOOST_CHECK_EQUAL(tree->GetEntries(), 3);
  rFile->Close();

  auto iFile =
      TFile::Open("RootMaterialMapIoGridGloballyIndexedTests.root", "READ");
  BOOST_REQUIRE(iFile != nullptr);
  // Default options - the persisted MaterialKind alone is enough to
  // reconstruct GloballyIndexed storage, no read-side configuration needed
  RootMaterialMapIo::Options readOptions;
  auto [surfaceMapsRead, volumeMapsRead] = accessor.read(*iFile, readOptions);
  BOOST_REQUIRE_EQUAL(surfaceMapsRead.size(), 2u);

  Vector2 lLow{0.5, 0.5};
  Vector2 lHigh{1.5, 0.5};

  const auto* readGrid0 = dynamic_cast<const GridSurfaceMaterial*>(
      surfaceMapsRead.at(geoID0).get());
  BOOST_REQUIRE(readGrid0 != nullptr);
  BOOST_CHECK(std::holds_alternative<GridSurfaceMaterial::GloballyIndexed>(
      readGrid0->storage()));
  BOOST_CHECK_CLOSE(readGrid0->materialSlab(lLow).material().X0(), 1., 1e-4);
  BOOST_CHECK_CLOSE(readGrid0->materialSlab(lHigh).material().X0(), 11., 1e-4);

  const auto* readGrid1 = dynamic_cast<const GridSurfaceMaterial*>(
      surfaceMapsRead.at(geoID1).get());
  BOOST_REQUIRE(readGrid1 != nullptr);
  BOOST_CHECK_CLOSE(readGrid1->materialSlab(lLow).material().X0(), 11., 1e-4);
  BOOST_CHECK_CLOSE(readGrid1->materialSlab(lHigh).material().X0(), 1., 1e-4);
}

// A second, distinct GloballyIndexed material vector in the same write
// session is a caller error - there is only ever one canonical global
// material vector per session, so its on-disk indices are always directly
// usable as absolute row numbers with no offset.
BOOST_AUTO_TEST_CASE(RootMaterialMapIoGridGloballyIndexedConflictThrows) {
  auto sharedMaterial0 = std::make_shared<std::vector<MaterialSlab>>();
  sharedMaterial0->emplace_back(Material::Vacuum(), 0.0);
  auto sharedMaterial1 = std::make_shared<std::vector<MaterialSlab>>();
  sharedMaterial1->emplace_back(Material::Vacuum(), 0.0);

  auto axisX = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 1., 1);
  auto axisY = IAxis::createEquidistant(AxisBoundaryType::Bound, 0., 1., 1);
  std::vector<std::vector<std::size_t>> indexPayload = {{0u}};

  auto gridMaterial0 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axisX, *axisY, sharedMaterial0, indexPayload);
  auto gridMaterial1 = GridSurfaceMaterialFactory::createGloballyIndexed(
      *axisX, *axisY, sharedMaterial1, indexPayload);

  auto geoID0 = GeometryIdentifier().withVolume(5).withSensitive(1);
  auto geoID1 = GeometryIdentifier().withVolume(5).withSensitive(2);

  auto rFile = TFile::Open("RootMaterialMapIoGridGloballyIndexedConflict.root",
                           "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  RootMaterialMapIo::Config cfg;
  RootMaterialMapIo accessor(cfg);
  RootMaterialMapIo::Options options;
  accessor.write(*rFile, geoID0, *gridMaterial0, options);
  BOOST_CHECK_THROW(accessor.write(*rFile, geoID1, *gridMaterial1, options),
                    std::invalid_argument);
  rFile->Close();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
