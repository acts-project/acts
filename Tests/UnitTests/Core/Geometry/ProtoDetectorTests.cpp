// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <sstream>

namespace {

/// @brief This method creates a world volume with
/// some sub structure, this detector is not yet
/// syncrhonized, it can then be typed into the
/// Detector or the TrackGeometry description
///
/// @return a proto world volume
Acts::ProtoDetector createProtoDetector() {
  // Container
  Acts::ProtoVolume detectorVolume;
  detectorVolume.name = "detector-container";
  detectorVolume.extent.set(Acts::binZ, -2000., 2000);

  // Beam Pipe volume
  Acts::ProtoVolume beamPipe;
  beamPipe.name = "beam-pipe";
  beamPipe.extent.set(Acts::binR, 0., 30.);

  // Pixel section
  Acts::ProtoVolume pixelContainer;
  pixelContainer.name = "pixel-container";
  pixelContainer.extent.set(Acts::binR, 40., 200);

  // Pixel volume sub structure
  Acts::ProtoVolume pixelNec;
  pixelNec.name = "pixel-nec";
  pixelNec.extent.set(Acts::binZ, -1900., -600);

  Acts::ProtoVolume pixelBarrel;
  pixelBarrel.name = "pixel-barrel";
  pixelBarrel.extent.set(Acts::binR, 41., 199.);
  pixelBarrel.extent.set(Acts::binZ, -550., 550.);

  Acts::ProtoVolume pixelBarrelL0;
  pixelBarrelL0.name = "pixel-barrel-l0";
  pixelBarrelL0.extent.set(Acts::binR, 45., 50.);
  pixelBarrelL0.layerType = Acts::Surface::SurfaceType::Cylinder;

  Acts::ProtoVolume pixelBarrelL1;
  pixelBarrelL1.name = "pixel-barrel-l1";
  pixelBarrelL1.extent.set(Acts::binR, 70., 80.);
  pixelBarrelL1.layerType = Acts::Surface::SurfaceType::Cylinder;

  pixelBarrel.constituentVolumes = {pixelBarrelL0, pixelBarrelL1};
  pixelBarrel.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1.})};

  Acts::ProtoVolume pixelPec;
  pixelPec.name = "pixel-pec";
  pixelPec.extent.set(Acts::binZ, 600., 1900.);

  pixelContainer.constituentVolumes = {pixelNec, pixelBarrel, pixelPec};
  pixelContainer.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};

  detectorVolume.constituentVolumes = {beamPipe, pixelContainer};
  detectorVolume.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1})};

  Acts::ProtoDetector detector;
  detector.name = "detector";
  detector.worldVolume = detectorVolume;

  return detector;
}

}  // namespace

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(ProtoTrackingGeometryTests) {
  // Get the raw proto detector description
  auto detector = createProtoDetector();
  detector.harmonize(true);

  // Get the top detector volume
  auto& detectorVolume = detector.worldVolume;

  // Check the container is NOT a layer Container
  BOOST_CHECK(detectorVolume.layerContainer == false);

  // The detector volume should have received maximum dimensions
  CHECK_CLOSE_ABS(detectorVolume.extent.min(Acts::binR), 0,
                  std::numeric_limits<ActsScalar>::epsilon());

  CHECK_CLOSE_ABS(detectorVolume.extent.max(Acts::binR), 200.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The detector cotainer should have binning in R
  BOOST_CHECK(detectorVolume.constituentBinning[0].type == Acts::arbitrary);
  BOOST_CHECK(detectorVolume.constituentBinning[0].binvalue == Acts::binR);

  const auto& binBoundaries = detectorVolume.constituentBinning[0].boundaries();
  BOOST_CHECK(binBoundaries.size() == 3u);
  CHECK_CLOSE_ABS(binBoundaries[0u], 0.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[1u], 35.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[2u], 200.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The first volume is the beam pipe, it should have gotten the
  // the z dimension
  auto& beamPipe = detectorVolume.constituentVolumes[0u];

  BOOST_CHECK(beamPipe.name == "beam-pipe");
  CHECK_CLOSE_ABS(beamPipe.extent.min(Acts::binZ), -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  CHECK_CLOSE_ABS(beamPipe.extent.max(Acts::binZ), 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The new beam pipe radius should have been applied
  CHECK_CLOSE_ABS(beamPipe.extent.max(Acts::binR), 35.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The second volume is the pixel detector
  auto& pixelContainer = detectorVolume.constituentVolumes[1u];
  BOOST_CHECK(pixelContainer.name == "pixel-container");
  BOOST_CHECK(pixelContainer.layerContainer == false);

  // Pixel contaienr should have fitting boundaries
  CHECK_CLOSE_ABS(pixelContainer.extent.min(Acts::binR), 35.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(Acts::binR), 200.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.min(Acts::binZ), -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(Acts::binZ), 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // All of the internal containers should now have synchronized
  // inner & outer boundaries
  for (auto& pv : pixelContainer.constituentVolumes) {
    CHECK_CLOSE_ABS(pv.extent.min(Acts::binR), 35.,
                    std::numeric_limits<ActsScalar>::epsilon());

    CHECK_CLOSE_ABS(pv.extent.max(Acts::binR), 200.,
                    std::numeric_limits<ActsScalar>::epsilon());
  }

  // The binning should have been estimated
  BOOST_CHECK(pixelContainer.constituentBinning[0].type == Acts::arbitrary);
  BOOST_CHECK(pixelContainer.constituentBinning[0].binvalue == Acts::binZ);

  const auto& binBoundariesZ =
      pixelContainer.constituentBinning[0].boundaries();
  BOOST_CHECK(binBoundariesZ.size() == 4u);
  CHECK_CLOSE_ABS(binBoundariesZ[0u], -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[1u], -575,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[2u], 575.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[3u], 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  auto& pixelNec = pixelContainer.constituentVolumes[0u];
  BOOST_CHECK(pixelNec.layerContainer == false);

  auto& pixelBarrel = pixelContainer.constituentVolumes[1u];
  BOOST_CHECK(pixelBarrel.layerContainer == true);

  auto& pixelPec = pixelContainer.constituentVolumes[2u];
  BOOST_CHECK(pixelPec.layerContainer == false);
}

BOOST_AUTO_TEST_CASE(ProtoDetectorTests) {
  // Get the raw proto detector description
  auto detector = createProtoDetector();
  detector.harmonize(false);

  std::cout << detector.toString() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
