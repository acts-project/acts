// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/Root/RootMaterialTrackPayload.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

using namespace Acts;

namespace {
// Make material surfaces
std::vector<std::shared_ptr<Acts::Surface>> materialSurfaces() {
  // Create surfaces
  std::vector<std::shared_ptr<Acts::Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 50.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 70.0,
                                           100.0)};
  // Assign geometry ids & material
  for (auto [is, surface] : enumerate(surfaces)) {
    // Assign unique geometry identifiers
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
    Acts::Material mat =
        Acts::Material::fromMolarDensity(100. + is, 33. + is, 3., 4., 5.);
    Acts::MaterialSlab mp(mat, 0.1);
    // Constructor from arguments
    auto msMaterial = std::make_shared<HomogeneousSurfaceMaterial>(mp, 1.);
    surface->assignSurfaceMaterial(msMaterial);
  }
  return surfaces;
}

// Make a material track
Acts::RecordedMaterialTrack materialTrack(
    const Acts::Vector3& position, const Acts::Vector3& momentum,
    const std::vector<std::shared_ptr<Acts::Surface>>& mSurfaces) {
  // Create a material track
  Acts::RecordedMaterialTrack rmTrack;
  // Fill the position and momentum
  rmTrack.first.first = position;
  rmTrack.first.second = momentum;
  // Fill the individual steps
  rmTrack.second.materialInX0 = 0.;
  rmTrack.second.materialInL0 = 0.;
  // Create a material interaction
  for (const auto& ms : mSurfaces) {
    // intersect the surface
    auto sfIntersection =
        ms->intersect(Acts::GeometryContext(), rmTrack.first.first,
                      rmTrack.first.second, Acts::BoundaryTolerance::None())
            .closest();
    // Create a material interaction
    Acts::MaterialInteraction mint;
    mint.surface = ms.get();
    mint.position = sfIntersection.position();
    mint.direction = rmTrack.first.second;
    mint.materialSlab =
        ms->surfaceMaterial()->materialSlab(sfIntersection.position());
    mint.pathCorrection = 1.0;
    rmTrack.second.materialInteractions.push_back(mint);
  }

  return rmTrack;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(RootPlugin)

BOOST_AUTO_TEST_CASE(RootMaterialTrackPayloadTests) {
  auto rFile =
      TFile::Open("RootMaterialTrackPayloadWriteTests.root", "RECREATE");

  // Some tree options

  // Create payload options
  auto rpTree =
      new TTree("material-tracks-plain", "Material Track Tree, Plain option");
  Acts::RootMaterialTrackPayload rmtpPlain;
  rmtpPlain.connectForWrite(*rpTree);

  auto rfTree =
      new TTree("material-tracks-full", "Material Track Tree, Plain option");
  Acts::RootMaterialTrackPayload rmtpFull(true, true, true, true, true);
  rmtpFull.connectForWrite(*rfTree);

  // Create the material surfaces
  auto mSurfaces = materialSurfaces();

  // Write 10 material tracks
  for (int i = 0; i < 10; i++) {
    auto position = Acts::Vector3(0.1 * i, 0.2 * i, 0.3 * i);
    auto momentum = Acts::Vector3(1. - 0.1 * i, 0.1 * i, 0.1 * i);
    auto rmTrack = materialTrack(position, momentum, mSurfaces);
    // Check the length
    BOOST_CHECK_EQUAL(rmTrack.second.materialInteractions.size(), 4u);
    // Fill the plain version
    rmtpPlain.write(GeometryContext(), i, rmTrack);
    rpTree->Fill();
    // Fill the full version
    rmtpFull.write(GeometryContext(), i, rmTrack);
    rfTree->Fill();
  }

  // Plain tree and full tree have 10 entries
  BOOST_CHECK_EQUAL(rpTree->GetEntries(), 10);
  BOOST_CHECK_EQUAL(rfTree->GetEntries(), 10);

  // Plain tree does not have pre/post step information
  BOOST_CHECK(rpTree->FindBranch("mat_sx") == nullptr);
  BOOST_CHECK(rpTree->FindBranch("mat_ex") == nullptr);

  BOOST_CHECK(rfTree->FindBranch("mat_sx") != nullptr);
  BOOST_CHECK(rfTree->FindBranch("mat_ex") != nullptr);

  rFile->Write();
  rFile->Close();

  // Read back in - plain infromation
  TChain rpChain("material-tracks-plain");
  rpChain.Add("RootMaterialTrackPayloadWriteTests.root");

  Acts::RootMaterialTrackPayload rmtpIn;
  rmtpIn.connectForRead(rpChain);

  for (int i = 0; i < 10; i++) {
    rpChain.GetEntry(i);
    auto rmTrack = rmtpIn.read();
    // All original are still here
    BOOST_CHECK_EQUAL(rmTrack.second.materialInteractions.size(), 4u);

    auto rPos = Acts::Vector3(0.1 * i, 0.2 * i, 0.3 * i);
    CHECK_CLOSE_ABS(rmTrack.first.first.x(), rPos.x(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.first.y(), rPos.y(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.first.z(), rPos.z(), 1e-6);

    auto rMom = Acts::Vector3(1. - 0.1 * i, 0.1 * i, 0.1 * i);
    CHECK_CLOSE_ABS(rmTrack.first.second.x(), rMom.x(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.second.y(), rMom.y(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.second.z(), rMom.z(), 1e-6);
  }

  // Read back in - full information
  TChain rfChain("material-tracks-full");
  rfChain.Add("RootMaterialTrackPayloadWriteTests.root");

  Acts::RootMaterialTrackPayload rmtpInFull;
  rmtpInFull.connectForRead(rfChain);

  for (int i = 0; i < 10; i++) {
    rfChain.GetEntry(i);
    auto rmTrack = rmtpInFull.read();
    // One is collapsed
    BOOST_CHECK_EQUAL(rmTrack.second.materialInteractions.size(), 3u);
    BOOST_CHECK(rmTrack.second.materialInX0 > 0.);

    auto rPos = Acts::Vector3(0.1 * i, 0.2 * i, 0.3 * i);
    CHECK_CLOSE_ABS(rmTrack.first.first.x(), rPos.x(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.first.y(), rPos.y(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.first.z(), rPos.z(), 1e-6);

    auto rMom = Acts::Vector3(1. - 0.1 * i, 0.1 * i, 0.1 * i);
    CHECK_CLOSE_ABS(rmTrack.first.second.x(), rMom.x(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.second.y(), rMom.y(), 1e-6);
    CHECK_CLOSE_ABS(rmTrack.first.second.z(), rMom.z(), 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
