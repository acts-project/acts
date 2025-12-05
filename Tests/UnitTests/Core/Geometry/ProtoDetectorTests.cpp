// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace Acts;

namespace ActsTests {

/// @brief This method creates a world volume with
/// some sub structure, this detector is not yet
/// synchronized, it can then be typed into the
/// Detector or the TrackGeometry description
///
/// @return a proto world volume
ProtoDetector createProtoDetector() {
  // Container
  ProtoVolume detectorVolume;
  detectorVolume.name = "detector-container";
  detectorVolume.extent.set(AxisDirection::AxisZ, -2000., 2000);

  // Beam Pipe volume
  ProtoVolume beamPipe;
  beamPipe.name = "beam-pipe";
  beamPipe.extent.set(AxisDirection::AxisR, 0., 30.);

  // Pixel section
  ProtoVolume pixelContainer;
  pixelContainer.name = "pixel-container";
  pixelContainer.extent.set(AxisDirection::AxisR, 40., 200);

  // Pixel volume sub structure
  ProtoVolume pixelNec;
  pixelNec.name = "pixel-nec";
  pixelNec.extent.set(AxisDirection::AxisZ, -1900., -600);

  ProtoVolume pixelBarrel;
  pixelBarrel.name = "pixel-barrel";
  pixelBarrel.extent.set(AxisDirection::AxisR, 41., 199.);
  pixelBarrel.extent.set(AxisDirection::AxisZ, -550., 550.);

  ProtoVolume pixelBarrelL0;
  pixelBarrelL0.name = "pixel-barrel-l0";
  pixelBarrelL0.extent.set(AxisDirection::AxisR, 45., 50.);
  pixelBarrelL0.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};

  ProtoVolume pixelBarrelL1;
  pixelBarrelL1.name = "pixel-barrel-l1";
  pixelBarrelL1.extent.set(AxisDirection::AxisR, 70., 80.);
  pixelBarrelL1.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};

  pixelBarrel.container = ProtoVolume::ContainerStructure{
      {pixelBarrelL0, pixelBarrelL1},
      {BinningData(open, AxisDirection::AxisR, {0., 1.})},
      true};

  ProtoVolume pixelPec;
  pixelPec.name = "pixel-pec";
  pixelPec.extent.set(AxisDirection::AxisZ, 600., 1900.);

  pixelContainer.container = ProtoVolume::ContainerStructure{
      {pixelNec, pixelBarrel, pixelPec},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})}};

  detectorVolume.container = ProtoVolume::ContainerStructure{
      {beamPipe, pixelContainer},
      {BinningData(open, AxisDirection::AxisR, {0., 1})}};

  ProtoDetector detector;
  detector.name = "detector";
  detector.worldVolume = detectorVolume;

  return detector;
}

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(ProtoTrackingGeometryTests) {
  // Get the raw proto detector description
  auto detector = createProtoDetector();
  detector.harmonize(true);

  // Get the top detector volume
  auto& detectorVolume = detector.worldVolume;

  // The detector volume should have received maximum dimensions
  CHECK_CLOSE_ABS(detectorVolume.extent.min(AxisDirection::AxisR), 0,
                  std::numeric_limits<double>::epsilon());

  CHECK_CLOSE_ABS(detectorVolume.extent.max(AxisDirection::AxisR), 200.,
                  std::numeric_limits<double>::epsilon());

  // The detector container should have binning in R
  BOOST_CHECK(detectorVolume.container.has_value());
  BOOST_CHECK(!detectorVolume.internal.has_value());

  auto& cts = detectorVolume.container.value();

  BOOST_CHECK_EQUAL(cts.constituentBinning.size(), 1u);
  BOOST_CHECK_EQUAL(cts.constituentBinning[0].type, arbitrary);
  BOOST_CHECK_EQUAL(cts.constituentBinning[0].binvalue, AxisDirection::AxisR);

  const auto& binBoundaries = cts.constituentBinning[0].boundaries();
  BOOST_CHECK_EQUAL(binBoundaries.size(), 3u);
  CHECK_CLOSE_ABS(binBoundaries[0u], 0.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[1u], 35.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(binBoundaries[2u], 200.,
                  std::numeric_limits<double>::epsilon());

  // The first volume is the beam pipe, it should have gotten the
  // z dimension
  auto& beamPipe = cts.constituentVolumes[0u];

  BOOST_CHECK_EQUAL(beamPipe.name, "beam-pipe");
  CHECK_CLOSE_ABS(beamPipe.extent.min(AxisDirection::AxisZ), -2000.,
                  std::numeric_limits<double>::epsilon());

  CHECK_CLOSE_ABS(beamPipe.extent.max(AxisDirection::AxisZ), 2000.,
                  std::numeric_limits<double>::epsilon());

  // The new beam pipe radius should have been applied
  CHECK_CLOSE_ABS(beamPipe.extent.max(AxisDirection::AxisR), 35.,
                  std::numeric_limits<double>::epsilon());

  // The second volume is the pixel detector
  auto& pixelContainer = cts.constituentVolumes[1u];
  BOOST_CHECK_EQUAL(pixelContainer.name, "pixel-container");

  // Pixel container should have fitting boundaries
  CHECK_CLOSE_ABS(pixelContainer.extent.min(AxisDirection::AxisR), 35.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(AxisDirection::AxisR), 200.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.min(AxisDirection::AxisZ), -2000.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(pixelContainer.extent.max(AxisDirection::AxisZ), 2000.,
                  std::numeric_limits<double>::epsilon());

  // The Pixel container has constituents
  BOOST_CHECK(pixelContainer.container.has_value());
  auto& cts1 = pixelContainer.container.value();

  // All of the internal containers should now have synchronized
  // inner & outer boundaries
  for (auto& pv : cts1.constituentVolumes) {
    CHECK_CLOSE_ABS(pv.extent.min(AxisDirection::AxisR), 35.,
                    std::numeric_limits<double>::epsilon());

    CHECK_CLOSE_ABS(pv.extent.max(AxisDirection::AxisR), 200.,
                    std::numeric_limits<double>::epsilon());
  }

  // The binning should have been estimated
  BOOST_CHECK_EQUAL(cts1.constituentBinning.size(), 1u);
  BOOST_CHECK_EQUAL(cts1.constituentBinning[0].type, arbitrary);
  BOOST_CHECK_EQUAL(cts1.constituentBinning[0].binvalue, AxisDirection::AxisZ);

  const auto& binBoundariesZ = cts1.constituentBinning[0].boundaries();
  BOOST_CHECK_EQUAL(binBoundariesZ.size(), 4u);
  CHECK_CLOSE_ABS(binBoundariesZ[0u], -2000.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[1u], -575,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[2u], 575.,
                  std::numeric_limits<double>::epsilon());
  CHECK_CLOSE_ABS(binBoundariesZ[3u], 2000.,
                  std::numeric_limits<double>::epsilon());

  // The second volume is the pixel barrel
  auto& pixelBarrel = cts1.constituentVolumes[1u];
  BOOST_CHECK_EQUAL(pixelBarrel.name, "pixel-barrel");

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

}  // namespace ActsTests
