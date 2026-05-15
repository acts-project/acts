// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsPlugins/Root/RootMeasurementIo.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <filesystem>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(RootMeasurementIoTestsWrite) {
  // Create temporary directory
  TemporaryDirectory tmpDir{};

  // Create the configuration
  std::vector<BoundIndices> recoIndices = {eBoundLoc0, eBoundLoc1};
  std::vector<BoundIndices> clusterIndices = {eBoundLoc0, eBoundLoc1};
  RootMeasurementIo::Config cfg{recoIndices, clusterIndices};
  RootMeasurementIo accessor(cfg);

  const std::filesystem::path filePath =
      tmpDir.path() / "RootMeasurementIoTests.root";
  // Create a test file
  auto rFile = TFile::Open(filePath.c_str(), "RECREATE");
  rFile->cd();
  BOOST_REQUIRE(rFile != nullptr);

  TTree measurementTree("measurements", "measurements");
  accessor.connectForWrite(measurementTree);

  // Fill some data to the geometry identifier
  auto geoId = GeometryIdentifier()
                   .withVolume(1)
                   .withLayer(2)
                   .withSensitive(4)
                   .withExtra(5);

  accessor.fillIdentification(1, geoId);
  accessor.fillTruthParameters(Vector2(0.1, 0.2), Vector4(1.1, 2.2, 3.3, 4.4),
                               Vector3(0.1, 0.2, 0.3),
                               std::make_pair(0.01, 0.02));
  accessor.fillBoundMeasurement({0.11, 0.22}, {0.01, 0.02}, {0, 1});
  accessor.fillGlobalPosition(Vector3(1.0, 2.0, 3.0));
  accessor.fillCluster(std::vector<std::tuple<int, int, float>>{
      {1, 2, 0.5}, {2, 3, 1.5}, {3, 4, 2.5}});

  measurementTree.Fill();
  measurementTree.Write();
  accessor.clear();
  rFile->Close();

  // Let's read it back
  rFile = TFile::Open(filePath.c_str(), "READ");
  BOOST_REQUIRE(rFile != nullptr);
  auto readTree = dynamic_cast<TTree*>(rFile->Get("measurements"));
  BOOST_REQUIRE(readTree != nullptr);

  BOOST_CHECK_EQUAL(readTree->GetEntries(), 1);
}

BOOST_AUTO_TEST_CASE(RootMeasurementIoExceptions) {
  // Create the configuration
  std::vector<BoundIndices> recoIndices = {eBoundLoc0};
  std::vector<BoundIndices> clusterIndices = {eBoundLoc0};
  RootMeasurementIo::Config cfg{recoIndices, clusterIndices};
  RootMeasurementIo accessor(cfg);

  BOOST_CHECK_THROW(accessor.fillBoundMeasurement({0.1, 0.2}, {0.01}, {0, 1}),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
