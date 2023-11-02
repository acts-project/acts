// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

/// @brief This method creates a world volume with
/// some sub structure, this detector is not yet
/// synchronized, it can then be typed into the
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
  pixelBarrelL0.internal = Acts::ProtoVolume::InternalStructure{
      Acts::Surface::SurfaceType::Cylinder};

  Acts::ProtoVolume pixelBarrelL1;
  pixelBarrelL1.name = "pixel-barrel-l1";
  pixelBarrelL1.extent.set(Acts::binR, 70., 80.);
  pixelBarrelL1.internal = Acts::ProtoVolume::InternalStructure{
      Acts::Surface::SurfaceType::Cylinder};

  pixelBarrel.container = Acts::ProtoVolume::ContainerStructure{
      {pixelBarrelL0, pixelBarrelL1},
      {Acts::BinningData(Acts::open, Acts::binR, {0., 1.})},
      true};

  Acts::ProtoVolume pixelPec;
  pixelPec.name = "pixel-pec";
  pixelPec.extent.set(Acts::binZ, 600., 1900.);

  pixelContainer.container = Acts::ProtoVolume::ContainerStructure{
      {pixelNec, pixelBarrel, pixelPec},
      {Acts::BinningData(Acts::open, Acts::binZ, {0., 1})}};

  detectorVolume.container = Acts::ProtoVolume::ContainerStructure{
      {beamPipe, pixelContainer},
      {Acts::BinningData(Acts::open, Acts::binR, {0., 1})}};

  Acts::ProtoDetector detector;
  detector.name = "detector";
  detector.worldVolume = detectorVolume;

  return detector;
}

}  // namespace

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(ProtoTrackingGeometryTests) {
  // Get the raw proto detector description
  auto detector = createProtoDetector();
  detector.harmonize(true);

  // Get the top detector volume
  auto& detectorVolume = detector.worldVolume;

  // The detector volume should have received maximum dimensions
  CHECK_CLOSE_ABS(detectorVolume.extent.min(Acts::binR), 0,
                  std::numeric_limits<ActsScalar>::epsilon());

  CHECK_CLOSE_ABS(detectorVolume.extent.max(Acts::binR), 200.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The detector container should have binning in R
  BOOST_CHECK(detectorVolume.container.has_value());
  BOOST_CHECK(not detectorVolume.internal.has_value());

  auto& cts = detectorVolume.container.value();

  BOOST_CHECK(cts.constituentBinning.size() == 1u);
  BOOST_CHECK(cts.constituentBinning[0].type == Acts::arbitrary);
  BOOST_CHECK(cts.constituentBinning[0].binvalue == Acts::binR);

  const auto& binBoundaries = cts.constituentBinning[0].boundaries();
  BOOST_CHECK(binBoundaries.size() == 3u);
  CHECK_CLOSE_ABS(binBoundaries[0u], 0.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[1u], 35.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[2u], 200.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The first volume is the beam pipe, it should have gotten the
  // z dimension
  auto& beamPipe = cts.constituentVolumes[0u];

  BOOST_CHECK(beamPipe.name == "beam-pipe");
  CHECK_CLOSE_ABS(beamPipe.extent.min(Acts::binZ), -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  CHECK_CLOSE_ABS(beamPipe.extent.max(Acts::binZ), 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The new beam pipe radius should have been applied
  CHECK_CLOSE_ABS(beamPipe.extent.max(Acts::binR), 35.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The second volume is the pixel detector
  auto& pixelContainer = cts.constituentVolumes[1u];
  BOOST_CHECK(pixelContainer.name == "pixel-container");

  // Pixel container should have fitting boundaries
  CHECK_CLOSE_ABS(pixelContainer.extent.min(Acts::binR), 35.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(Acts::binR), 200.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.min(Acts::binZ), -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(Acts::binZ), 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The Pixel container has constituents
  BOOST_CHECK(pixelContainer.container.has_value());
  auto& cts1 = pixelContainer.container.value();

  // All of the internal containers should now have synchronized
  // inner & outer boundaries
  for (auto& pv : cts1.constituentVolumes) {
    CHECK_CLOSE_ABS(pv.extent.min(Acts::binR), 35.,
                    std::numeric_limits<ActsScalar>::epsilon());

    CHECK_CLOSE_ABS(pv.extent.max(Acts::binR), 200.,
                    std::numeric_limits<ActsScalar>::epsilon());
  }

  // The binning should have been estimated
  BOOST_CHECK(cts1.constituentBinning.size() == 1u);
  BOOST_CHECK(cts1.constituentBinning[0].type == Acts::arbitrary);
  BOOST_CHECK(cts1.constituentBinning[0].binvalue == Acts::binZ);

  const auto& binBoundariesZ = cts1.constituentBinning[0].boundaries();
  BOOST_CHECK(binBoundariesZ.size() == 4u);
  CHECK_CLOSE_ABS(binBoundariesZ[0u], -2000.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[1u], -575,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[2u], 575.,
                  std::numeric_limits<ActsScalar>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[3u], 2000.,
                  std::numeric_limits<ActsScalar>::epsilon());

  // The second volume is the pixel barrel
  auto& pixelBarrel = cts1.constituentVolumes[1u];
  BOOST_CHECK(pixelBarrel.name == "pixel-barrel");

  // It is a container volume value
  BOOST_CHECK(pixelBarrel.container.has_value());
  auto& cts2 = pixelBarrel.container.value();
  // It is, however, a layer container
  BOOST_CHECK(cts2.layerContainer);
  for (auto& lVolume : cts2.constituentVolumes) {
    BOOST_CHECK(lVolume.internal.has_value());
  }
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
