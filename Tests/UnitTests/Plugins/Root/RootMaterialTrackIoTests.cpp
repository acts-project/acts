// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Plugins/Root/RootMaterialTrackIo.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <memory>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

using namespace Acts;

std::vector<std::shared_ptr<Surface>> createTestSurfaces() {
  std::vector<std::shared_ptr<Surface>> surfaces;
  Transform3 transform = Transform3::Identity();
  // Create some dummy surfaces for testing
  for (int i = 0; i < 10; ++i) {
    // Create a surface
    transform.translation() = Vector3(0., 0., i * 100.);
    auto surface = Surface::makeShared<PlaneSurface>(transform);
    // Assign some material
    // construct the material properties from arguments
    Material mat = Material::fromMolarDensity(
        1. + i * 0.5, 2. + i * 0.5, 3. + i * 0.5, 4. + i * 0.5, 5. + i * 0.5);
    MaterialSlab mp(mat, 0.1);
    auto hMaterial = std::make_shared<HomogeneousSurfaceMaterial>(mp);
    surface->assignSurfaceMaterial(hMaterial);
    // Assign the geometry identifier
    GeometryIdentifier id =
        GeometryIdentifier().withVolume(i * 1).withSensitive(i);
    surface->assignGeometryId(id);
    surfaces.push_back(surface);
  }
  return surfaces;
}

BOOST_AUTO_TEST_SUITE(RootMaterialTrackIoTests)

BOOST_AUTO_TEST_CASE(RootReadWriteMaterialTracks) {
  // Create the setup first
  auto surfaces = createTestSurfaces();
  auto tContext = GeometryContext();
  // Fake some material tracks
  std::vector<RecordedMaterialTrack> rmTracks;
  Vector3 position(0.0, 0.0, 0.0);
  for (std::size_t it = 0; it < 10; ++it) {
    RecordedMaterial rMaterial;
    Vector3 direction = Vector3(0.005 * it, 0.005 * it, 1.0).normalized();
    for (const auto& surface : surfaces) {
      // Create a material interaction
      auto sIntersection = surface->intersect(tContext, position, direction);
      for (auto& si : sIntersection.split()) {
        if (si.isValid()) {
          MaterialInteraction interaction;
          interaction.position = si.position();
          interaction.direction = direction;
          interaction.intersection = si.position();
          interaction.intersectionID = surface->geometryId();
          interaction.time = 0.0;  // Set a dummy time
          interaction.surface = surface.get();
          interaction.materialSlab =
              surface->surfaceMaterial()->materialSlab(position);
          interaction.pathCorrection = 1. / std::abs(direction.z());
          rMaterial.materialInteractions.emplace_back(interaction);
          rMaterial.materialInX0 += interaction.materialSlab.thicknessInX0() *
                                    interaction.pathCorrection;
          rMaterial.materialInL0 += interaction.materialSlab.thicknessInL0() *
                                    interaction.pathCorrection;
        }
      }
    }
    RecordedMaterialTrack rmTrack{{position, direction}, rMaterial};
    rmTracks.push_back(rmTrack);
  }

  // Create the IO object - standard writing
  RootMaterialTrackIo::Config rmtConfig;
  RootMaterialTrackIo rmtIO(rmtConfig);

  auto materialFile = TFile::Open("MaterialTracksDefault.root", "RECREATE");
  TTree materialTree("materialTree", "Material Track Tree");
  rmtIO.connectForWrite(materialTree);

  // Write the material tracks
  unsigned int imt = 0;
  for (const auto& rmTrack : rmTracks) {
    rmtIO.write(tContext, imt++, rmTrack);
    materialTree.Fill();
  }
  materialTree.Write();
  materialFile->Close();

  // Read the material tracks back
  auto materialChain = TChain("materialTree");
  materialChain.Add("MaterialTracksDefault.root");
  RootMaterialTrackIo rmtIORead(rmtConfig);
  rmtIORead.connectForRead(materialChain);

  std::vector<RecordedMaterialTrack> readTracks;
  for (unsigned int i = 0; i < materialChain.GetEntries(); ++i) {
    materialChain.GetEntry(i);
    auto rmTrackRead = rmtIORead.read();
    readTracks.push_back(rmTrackRead);
  }
  BOOST_REQUIRE_EQUAL(rmTracks.size(), readTracks.size());
  for (std::size_t i = 0; i < rmTracks.size(); ++i) {
    const auto& rmTrack = rmTracks[i];
    const auto& rmTrackRead = readTracks[i];
    const auto& rPos = rmTrack.first.first;
    const auto& rDir = rmTrack.first.second;
    const auto& rPosRead = rmTrackRead.first.first;
    const auto& rDirRead = rmTrackRead.first.second;

    CHECK_CLOSE_ABS(rPosRead.x(), rPos.x(), 1e-6);
    CHECK_CLOSE_ABS(rPosRead.y(), rPos.y(), 1e-6);
    CHECK_CLOSE_ABS(rPosRead.z(), rPos.z(), 1e-6);
    CHECK_CLOSE_ABS(rDirRead.x(), rDir.x(), 1e-6);
    CHECK_CLOSE_ABS(rDirRead.y(), rDir.y(), 1e-6);
    CHECK_CLOSE_ABS(rDirRead.z(), rDir.z(), 1e-6);

    CHECK_CLOSE_ABS(rmTrackRead.second.materialInX0,
                    rmTrack.second.materialInX0, 1e-6);
    CHECK_CLOSE_ABS(rmTrackRead.second.materialInL0,
                    rmTrack.second.materialInL0, 1e-6);
    BOOST_REQUIRE_EQUAL(rmTrackRead.second.materialInteractions.size(),
                        rmTrack.second.materialInteractions.size());
    for (std::size_t j = 0; j < rmTrackRead.second.materialInteractions.size();
         ++j) {
      const auto& interaction = rmTrack.second.materialInteractions[j];
      const auto& interactionRead = rmTrackRead.second.materialInteractions[j];
      CHECK_CLOSE_ABS(interactionRead.position.x(), interaction.position.x(),
                      1e-6);
      CHECK_CLOSE_ABS(interactionRead.position.y(), interaction.position.y(),
                      1e-6);
      CHECK_CLOSE_ABS(interactionRead.position.z(), interaction.position.z(),
                      1e-6);
      CHECK_CLOSE_ABS(interactionRead.direction.x(), interaction.direction.x(),
                      1e-6);
      CHECK_CLOSE_ABS(interactionRead.direction.y(), interaction.direction.y(),
                      1e-6);
      CHECK_CLOSE_ABS(interactionRead.direction.z(), interaction.direction.z(),
                      1e-6);
    }
  }

  // Create the IO object - full writing
  RootMaterialTrackIo::Config rmtConfigFull;
  rmtConfigFull.prePostStepInfo = true;
  rmtConfigFull.surfaceInfo = true;
  rmtConfigFull.volumeInfo = true;
  RootMaterialTrackIo rmtIOFull(rmtConfigFull);

  auto materialFileFull = TFile::Open("MaterialTracksFull.root", "RECREATE");
  TTree materialTreeFull("materialTree", "Material Track Tree");
  rmtIOFull.connectForWrite(materialTreeFull);

  // Write the material tracks
  imt = 0;
  for (const auto& rmTrack : rmTracks) {
    rmtIOFull.write(tContext, imt++, rmTrack);
    materialTreeFull.Fill();
  }
  materialTreeFull.Write();
  materialFileFull->Close();

  // Read the material tracks back
  auto materialChainFull = TChain("materialTree");
  materialChainFull.Add("MaterialTracksFull.root");
  RootMaterialTrackIo rmtIOReadFull(rmtConfigFull);
  rmtIOReadFull.connectForRead(materialChainFull);

  std::vector<RecordedMaterialTrack> readTracksFull;
  for (unsigned int i = 0; i < materialChainFull.GetEntries(); ++i) {
    materialChainFull.GetEntry(i);
    auto rmTrackRead = rmtIOReadFull.read();
    readTracksFull.push_back(rmTrackRead);
  }
  BOOST_REQUIRE_EQUAL(rmTracks.size(), readTracksFull.size());

  for (std::size_t i = 0; i < rmTracks.size(); ++i) {
    const auto& rmTrack = rmTracks[i];
    const auto& rmTrackRead = readTracksFull[i];

    CHECK_CLOSE_ABS(rmTrackRead.second.materialInX0,
                    rmTrack.second.materialInX0, 1e-6);
    CHECK_CLOSE_ABS(rmTrackRead.second.materialInL0,
                    rmTrack.second.materialInL0, 1e-6);
    BOOST_REQUIRE_EQUAL(rmTrackRead.second.materialInteractions.size(),
                        rmTrack.second.materialInteractions.size());

    for (std::size_t j = 0; j < rmTrackRead.second.materialInteractions.size();
         ++j) {
      const auto& interaction = rmTrack.second.materialInteractions[j];
      const auto& interactionRead = rmTrackRead.second.materialInteractions[j];

      CHECK_CLOSE_ABS(interactionRead.intersection.x(),
                      interaction.intersection.x(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.intersection.y(),
                      interaction.intersection.y(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.intersection.z(),
                      interaction.intersection.z(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.time, interaction.time, 1e-6);
      CHECK_CLOSE_ABS(interactionRead.pathCorrection,
                      interaction.pathCorrection, 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.thicknessInX0(),
                      interaction.materialSlab.thicknessInX0(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.thicknessInL0(),
                      interaction.materialSlab.thicknessInL0(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.material().X0(),
                      interaction.materialSlab.material().X0(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.material().L0(),
                      interaction.materialSlab.material().L0(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.material().Ar(),
                      interaction.materialSlab.material().Ar(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.material().Z(),
                      interaction.materialSlab.material().Z(), 1e-6);
      CHECK_CLOSE_ABS(interactionRead.materialSlab.material().massDensity(),
                      interaction.materialSlab.material().massDensity(), 1e-6);
      BOOST_CHECK_EQUAL(interactionRead.intersectionID.value(),
                        interaction.intersectionID.value());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
