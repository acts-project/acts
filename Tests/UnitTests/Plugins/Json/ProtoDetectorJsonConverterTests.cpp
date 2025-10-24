// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsPlugins/Json/ProtoDetectorJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iterator>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "EqualityHelpers.hpp"

using namespace Acts;

namespace ActsTests {

/// @brief Helper method to compare proto volumes
/// @param one the first volume object
/// @param two the second volume object
/// @param tolerance the tolerance
/// @return a boolean to see if they are equal
bool isEqual(const ProtoVolume& one, const ProtoVolume& two,
             const double tolerance = 0.) {
  bool nameEq = (one.name == two.name);
  // Name
  BOOST_CHECK(nameEq);
  // Extent
  bool extentEq = isEqual(one.extent, two.extent, tolerance);
  BOOST_CHECK(extentEq);

  // Check internal structure
  bool internalValueEq = (one.internal.has_value() == two.internal.has_value());
  BOOST_CHECK(internalValueEq);
  bool internalEq = internalValueEq;
  if (one.internal.has_value() && two.internal.has_value()) {
    // Check consistency of the internal structure
    const auto& itsOne = one.internal.value();
    const auto& itsTwo = two.internal.value();
    bool layerTypeEq = (itsOne.layerType == itsTwo.layerType);
    BOOST_CHECK(layerTypeEq);
    internalEq = layerTypeEq;
    bool sBinningSizeEq =
        (itsOne.surfaceBinning.size() == itsTwo.surfaceBinning.size());
    BOOST_CHECK(sBinningSizeEq);
    internalEq = internalEq && sBinningSizeEq;
    for (auto [isb, sb] : enumerate(itsOne.surfaceBinning)) {
      bool sBinningEq = isEqual(sb, itsTwo.surfaceBinning[isb], tolerance);
      BOOST_CHECK(sBinningEq);
      internalEq = internalEq && sBinningEq;
    }
  }
  BOOST_CHECK(internalEq);

  // Check container structure
  bool containerValueEq =
      (one.container.has_value() == two.container.has_value());
  BOOST_CHECK(containerValueEq);
  bool containerEq = containerValueEq;
  if (one.container.has_value() && two.container.has_value()) {
    // Check consistency of the container structure
    const auto& ctsOne = one.container.value();
    const auto& ctsTwo = two.container.value();
    bool layerContainerEq = (ctsOne.layerContainer == ctsTwo.layerContainer);
    BOOST_CHECK(layerContainerEq);
    containerEq = layerContainerEq;
    bool cBinningSizeEq =
        ctsOne.constituentBinning.size() == ctsTwo.constituentBinning.size();
    containerEq = containerEq && cBinningSizeEq;
    BOOST_CHECK(cBinningSizeEq);
    for (auto [icb, cb] : enumerate(ctsOne.constituentBinning)) {
      bool cBinningEq = isEqual(cb, ctsTwo.constituentBinning[icb], tolerance);
      BOOST_CHECK(cBinningEq);
      containerEq = containerEq && cBinningEq;
    }
    // Recursively walk down
    bool cSizeEq =
        (ctsOne.constituentVolumes.size() == ctsTwo.constituentVolumes.size());
    BOOST_CHECK(cSizeEq);
    containerEq = cSizeEq;
    for (auto [ic, cOne] : enumerate(ctsOne.constituentVolumes)) {
      const auto& cTwo = ctsTwo.constituentVolumes[ic];
      bool cEq = isEqual(cOne, cTwo, tolerance);
      BOOST_CHECK(cEq);
      containerEq = containerEq && cEq;
    }
  }
  BOOST_CHECK(containerEq);

  // Give the overall judgement
  return nameEq && extentEq && internalEq && containerEq;
}

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(ProtoDetectorRoundTrip) {
  ExtentEnvelope cylinderLayerEnvelope = ExtentEnvelope::Zero();
  cylinderLayerEnvelope[AxisDirection::AxisR] = {1., 1.};
  cylinderLayerEnvelope[AxisDirection::AxisZ] = {2., 2.};

  ExtentEnvelope discLayerEnvelope = ExtentEnvelope::Zero();
  discLayerEnvelope[AxisDirection::AxisR] = {1., 1.};
  discLayerEnvelope[AxisDirection::AxisZ] = {1., 1.};

  // Beam Pipe container
  ProtoVolume beamPipeContainer;
  beamPipeContainer.name = "odd-beam-pipe";
  beamPipeContainer.extent.set(AxisDirection::AxisR, 0., 25);
  ProtoVolume beamPipe;
  beamPipe.name = "odd-beam-pipe-l";
  beamPipe.extent.set(AxisDirection::AxisR, 2., 24.);
  beamPipe.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  beamPipeContainer.container = ProtoVolume::ContainerStructure{
      {beamPipe}, {BinningData(open, AxisDirection::AxisR, {0., 1.})}, true};

  // Pixel section
  ProtoVolume pixelContainer;
  pixelContainer.name = "odd-pixel";
  pixelContainer.extent.set(AxisDirection::AxisR, 25., 200);

  ProtoVolume pixelNec;
  pixelNec.name = "odd-pixel-nec";
  pixelNec.extent.set(AxisDirection::AxisZ, -3100., -580);

  ProtoVolume pixNecD6;
  pixNecD6.name = "odd-pixel-nec-d6";
  pixNecD6.extent.set(AxisDirection::AxisZ, -1540., -1500);
  ProtoVolume pixNecD5;
  pixNecD5.name = "odd-pixel-nec-d5";
  pixNecD5.extent.set(AxisDirection::AxisZ, -1340., -1300);
  ProtoVolume pixNecD4;
  pixNecD4.name = "odd-pixel-nec-d4";
  pixNecD4.extent.set(AxisDirection::AxisZ, -1140., -1100);
  ProtoVolume pixNecD3;
  pixNecD3.name = "odd-pixel-nec-d3";
  pixNecD3.extent.set(AxisDirection::AxisZ, -1000., -960.);
  ProtoVolume pixNecD2;
  pixNecD2.name = "odd-pixel-nec-d2";
  pixNecD2.extent.set(AxisDirection::AxisZ, -860., -820);
  ProtoVolume pixNecD1;
  pixNecD1.name = "odd-pixel-nec-d1";
  pixNecD1.extent.set(AxisDirection::AxisZ, -740., -700);
  ProtoVolume pixNecD0;
  pixNecD0.name = "odd-pixel-nec-d0";
  pixNecD0.extent.set(AxisDirection::AxisZ, -640., -600);
  pixelNec.container = ProtoVolume::ContainerStructure{
      {pixNecD6, pixNecD5, pixNecD4, pixNecD3, pixNecD2, pixNecD1, pixNecD0},
      {BinningData(open, AxisDirection::AxisZ, {0., 1.})},
      true};

  BinningData pixEcBinningR =
      BinningData(open, AxisDirection::AxisR, 2., 0., 1.);
  BinningData pixEcBinningPhi = BinningData(
      closed, AxisDirection::AxisPhi, 30., -std::numbers::pi, std::numbers::pi);

  for (auto& cv : pixelNec.container.value().constituentVolumes) {
    cv.extent.setEnvelope(discLayerEnvelope);
    cv.internal = ProtoVolume::InternalStructure{
        Surface::SurfaceType::Disc, {pixEcBinningR, pixEcBinningPhi}};
  }

  ProtoVolume pixelBarrel;
  pixelBarrel.name = "odd-pixel-barrel";
  pixelBarrel.extent.set(AxisDirection::AxisZ, -580., 580);

  ProtoVolume pixBarrelL0;
  pixBarrelL0.name = "odd-pixel-barrel-l0";
  pixBarrelL0.extent.set(AxisDirection::AxisR, 28., 48.);
  pixBarrelL0.extent.set(AxisDirection::AxisZ, -580., 580);
  ProtoVolume pixBarrelL1;
  pixBarrelL1.name = "odd-pixel-barrel-l1";
  pixBarrelL1.extent.set(AxisDirection::AxisR, 62., 76);
  pixBarrelL1.extent.set(AxisDirection::AxisZ, -580., 580);
  ProtoVolume pixBarrelL2;
  pixBarrelL2.name = "odd-pixel-barrel-l2";
  pixBarrelL2.extent.set(AxisDirection::AxisR, 100., 120.);
  pixBarrelL2.extent.set(AxisDirection::AxisZ, -580., 580);
  ProtoVolume pixBarrelL3;
  pixBarrelL3.name = "odd-pixel-barrel-l3";
  pixBarrelL3.extent.set(AxisDirection::AxisR, 160., 180.);
  pixBarrelL3.extent.set(AxisDirection::AxisZ, -580., 580);

  pixelBarrel.container = ProtoVolume::ContainerStructure{
      {pixBarrelL0, pixBarrelL1, pixBarrelL2, pixBarrelL3},
      {BinningData(open, AxisDirection::AxisR, {0., 1})},
      true};

  for (auto& cv : pixelBarrel.container.value().constituentVolumes) {
    cv.extent.setEnvelope(cylinderLayerEnvelope);
    cv.internal =
        ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  }

  ProtoVolume pixelPec;
  pixelPec.name = "odd-pixel-pec";
  pixelPec.extent.set(AxisDirection::AxisZ, 580., 3100.);

  ProtoVolume pixPecD0;
  pixPecD0.name = "odd-pixel-pec-d0";
  pixPecD0.extent.set(AxisDirection::AxisZ, 600., 640);
  ProtoVolume pixPecD1;
  pixPecD1.name = "odd-pixel-pec-d1";
  pixPecD1.extent.set(AxisDirection::AxisZ, 700., 740);
  ProtoVolume pixPecD2;
  pixPecD2.name = "odd-pixel-pec-d2";
  pixPecD2.extent.set(AxisDirection::AxisZ, 820., 860.);
  ProtoVolume pixPecD3;
  pixPecD3.name = "odd-pixel-pec-d3";
  pixPecD3.extent.set(AxisDirection::AxisZ, 960., 1000.);
  ProtoVolume pixPecD4;
  pixPecD4.name = "odd-pixel-pec-d4";
  pixPecD4.extent.set(AxisDirection::AxisZ, 1100., 1140);
  ProtoVolume pixPecD5;
  pixPecD5.name = "odd-pixel-pec-d5";
  pixPecD5.extent.set(AxisDirection::AxisZ, 1300., 1340.);
  ProtoVolume pixPecD6;
  pixPecD6.name = "odd-pixel-pec-d6";
  pixPecD6.extent.set(AxisDirection::AxisZ, 1500., 1540.);

  pixelPec.container = ProtoVolume::ContainerStructure{
      {pixPecD0, pixPecD1, pixPecD2, pixPecD3, pixPecD4, pixPecD5, pixPecD6},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};

  for (auto& cv : pixelPec.container.value().constituentVolumes) {
    cv.extent.setEnvelope(discLayerEnvelope);
    cv.internal = ProtoVolume::InternalStructure{
        Surface::SurfaceType::Disc, {pixEcBinningR, pixEcBinningPhi}};
  }

  pixelContainer.container = ProtoVolume::ContainerStructure{
      {pixelNec, pixelBarrel, pixelPec},
      {BinningData(open, AxisDirection::AxisZ, {-3100., -580., 580., 3100.})}};

  // Short Strip section
  ProtoVolume pstContainer;
  pstContainer.name = "odd-pst";
  pstContainer.extent.set(AxisDirection::AxisR, 200., 210.);
  ProtoVolume pst;
  pst.name = "odd-pst-l";
  pst.extent.set(AxisDirection::AxisR, 201., 209.);
  pst.internal = ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  pstContainer.container = ProtoVolume::ContainerStructure{
      {pst}, {BinningData(open, AxisDirection::AxisR, {0., 1.})}, true};

  // Short Strip section
  ProtoVolume sstripContainer;
  sstripContainer.name = "odd-sstrip";
  sstripContainer.extent.set(AxisDirection::AxisR, 210., 720);

  BinningData sstripEcBinningR =
      BinningData(open, AxisDirection::AxisR, 3., 0., 1.);
  BinningData sstripEcBinningPhi = BinningData(
      closed, AxisDirection::AxisPhi, 42., -std::numbers::pi, std::numbers::pi);

  ProtoVolume sstripNec;
  sstripNec.name = "odd-sstrip-nec";
  sstripNec.extent.set(AxisDirection::AxisZ, -3100., -1200);
  ProtoVolume sstripNecD5;
  sstripNecD5.name = "odd-sstrip-nec-d5";
  sstripNecD5.extent.set(AxisDirection::AxisZ, -3000, -2900.);
  ProtoVolume sstripNecD4;
  sstripNecD4.name = "odd-sstrip-nec-d4";
  sstripNecD4.extent.set(AxisDirection::AxisZ, -2600., -2500.);
  ProtoVolume sstripNecD3;
  sstripNecD3.name = "odd-sstrip-nec-d3";
  sstripNecD3.extent.set(AxisDirection::AxisZ, -2250, -2150.);
  ProtoVolume sstripNecD2;
  sstripNecD2.name = "odd-sstrip-nec-d2";
  sstripNecD2.extent.set(AxisDirection::AxisZ, -1900, -1800.);
  ProtoVolume sstripNecD1;
  sstripNecD1.name = "odd-sstrip-nec-d1";
  sstripNecD1.extent.set(AxisDirection::AxisZ, -1600., -1500.);
  ProtoVolume sstripNecD0;
  sstripNecD0.name = "odd-sstrip-nec-d0";
  sstripNecD0.extent.set(AxisDirection::AxisZ, -1350., -1250.);

  sstripNec.container = ProtoVolume::ContainerStructure{
      {sstripNecD5, sstripNecD4, sstripNecD3, sstripNecD2, sstripNecD1,
       sstripNecD0},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};

  for (auto& cv : sstripNec.container.value().constituentVolumes) {
    cv.extent.setEnvelope(discLayerEnvelope);
    cv.internal = ProtoVolume::InternalStructure{
        Surface::SurfaceType::Disc, {sstripEcBinningR, sstripEcBinningPhi}};
  }

  ProtoVolume sstripBarrel;
  sstripBarrel.name = "odd-sstrip-barrel";
  sstripBarrel.extent.set(AxisDirection::AxisZ, -1200., 1200);

  ProtoVolume sstripBarrelL0;
  sstripBarrelL0.name = "odd-sstrip-barrel-l0";
  sstripBarrelL0.extent.set(AxisDirection::AxisR, 240., 280.);
  ProtoVolume sstripBarrelL1;
  sstripBarrelL1.name = "odd-sstrip-barrel-l1";
  sstripBarrelL1.extent.set(AxisDirection::AxisR, 340., 380.);
  ProtoVolume sstripBarrelL2;
  sstripBarrelL2.name = "odd-sstrip-barrel-l2";
  sstripBarrelL2.extent.set(AxisDirection::AxisR, 480., 520.);
  ProtoVolume sstripBarrelL3;
  sstripBarrelL3.name = "odd-sstrip-barrel-l3";
  sstripBarrelL3.extent.set(AxisDirection::AxisR, 640., 680.);

  sstripBarrel.container = ProtoVolume::ContainerStructure{
      {sstripBarrelL0, sstripBarrelL1, sstripBarrelL2, sstripBarrelL3},
      {BinningData(open, AxisDirection::AxisR, {0., 1})},
      true};

  for (auto& cv : sstripBarrel.container.value().constituentVolumes) {
    cv.extent.setEnvelope(cylinderLayerEnvelope);
    cv.internal =
        ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  }

  ProtoVolume sstripPec;
  sstripPec.name = "odd-sstrip-pec";
  sstripPec.extent.set(AxisDirection::AxisZ, 1200., 3100);

  ProtoVolume sstripPecD0;
  sstripPecD0.name = "odd-sstrip-pec-d0";
  sstripPecD0.extent.set(AxisDirection::AxisZ, 1250., 1350);
  ProtoVolume sstripPecD1;
  sstripPecD1.name = "odd-sstrip-pec-d1";
  sstripPecD1.extent.set(AxisDirection::AxisZ, 1500., 1600.);
  ProtoVolume sstripPecD2;
  sstripPecD2.name = "odd-sstrip-pec-d2";
  sstripPecD2.extent.set(AxisDirection::AxisZ, 1800., 1900.);
  ProtoVolume sstripPecD3;
  sstripPecD3.name = "odd-sstrip-pec-d3";
  sstripPecD3.extent.set(AxisDirection::AxisZ, 2150., 2250.);
  ProtoVolume sstripPecD4;
  sstripPecD4.name = "odd-sstrip-pec-d4";
  sstripPecD4.extent.set(AxisDirection::AxisZ, 2500., 2600.);
  ProtoVolume sstripPecD5;
  sstripPecD5.name = "odd-sstrip-pec-d5";
  sstripPecD5.extent.set(AxisDirection::AxisZ, 2900., 3000.);

  sstripPec.container = ProtoVolume::ContainerStructure{
      {sstripPecD0, sstripPecD1, sstripPecD2, sstripPecD3, sstripPecD4,
       sstripPecD5},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};
  for (auto& cv : sstripPec.container.value().constituentVolumes) {
    cv.extent.setEnvelope(discLayerEnvelope);
    cv.internal = ProtoVolume::InternalStructure{
        Surface::SurfaceType::Disc, {sstripEcBinningR, sstripEcBinningPhi}};
  }

  sstripContainer.container = ProtoVolume::ContainerStructure{
      {sstripNec, sstripBarrel, sstripPec},
      {BinningData(open, AxisDirection::AxisZ,
                   {-3100., -1200., 1200., 3100.})}};

  // Long Strip section
  ProtoVolume lstripContainer;
  lstripContainer.name = "odd-lstrip";
  lstripContainer.extent.set(AxisDirection::AxisR, 720, 1100.);

  ProtoVolume lstripNec;
  lstripNec.name = "odd-lstrip-nec";
  lstripNec.extent.set(AxisDirection::AxisZ, -3100., -1200);
  ProtoVolume lstripNecD5;
  lstripNecD5.name = "odd-lstrip-nec-d5";
  lstripNecD5.extent.set(AxisDirection::AxisZ, -3050, -2900.);
  ProtoVolume lstripNecD4;
  lstripNecD4.name = "odd-lstrip-nec-d4";
  lstripNecD4.extent.set(AxisDirection::AxisZ, -2650., -2500.);
  ProtoVolume lstripNecD3;
  lstripNecD3.name = "odd-lstrip-nec-d3";
  lstripNecD3.extent.set(AxisDirection::AxisZ, -2300, -2150.);
  ProtoVolume lstripNecD2;
  lstripNecD2.name = "odd-lstrip-nec-d2";
  lstripNecD2.extent.set(AxisDirection::AxisZ, -1950, -1800.);
  ProtoVolume lstripNecD1;
  lstripNecD1.name = "odd-lstrip-nec-d1";
  lstripNecD1.extent.set(AxisDirection::AxisZ, -1650., -1500.);
  ProtoVolume lstripNecD0;
  lstripNecD0.name = "odd-lstrip-nec-d0";
  lstripNecD0.extent.set(AxisDirection::AxisZ, -1400., -1250.);

  lstripNec.container = ProtoVolume::ContainerStructure{
      {lstripNecD5, lstripNecD4, lstripNecD3, lstripNecD2, lstripNecD1,
       lstripNecD0},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};

  for (auto& cv : lstripNec.container.value().constituentVolumes) {
    cv.extent.setEnvelope(discLayerEnvelope);
    cv.internal = ProtoVolume::InternalStructure{Surface::SurfaceType::Disc};
  }

  ProtoVolume lstripBarrel;
  lstripBarrel.name = "odd-lstrip-barrel";
  lstripBarrel.extent.set(AxisDirection::AxisZ, -1200., 1200);

  ProtoVolume lstripBarrelL0;
  lstripBarrelL0.name = "odd-lstrip-barrel-l0";
  lstripBarrelL0.extent.set(AxisDirection::AxisR, 800., 840.);
  ProtoVolume lstripBarrelL1;
  lstripBarrelL1.name = "odd-lstrip-barrel-l1";
  lstripBarrelL1.extent.set(AxisDirection::AxisR, 1000., 1050.);

  lstripBarrel.container = ProtoVolume::ContainerStructure{
      {lstripBarrelL0, lstripBarrelL1},
      {BinningData(open, AxisDirection::AxisR, {0., 1})},
      true};

  for (auto& cv : lstripBarrel.container.value().constituentVolumes) {
    cv.extent.setEnvelope(cylinderLayerEnvelope);
    cv.internal =
        ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  }

  ProtoVolume lstripPec;
  lstripPec.name = "odd-lstrip-pec";
  lstripPec.extent.set(AxisDirection::AxisZ, 1200., 3100);

  ProtoVolume lstripPecD0;
  lstripPecD0.name = "odd-lstrip-pec-d0";
  lstripPecD0.extent.set(AxisDirection::AxisZ, 1250., 1400);
  ProtoVolume lstripPecD1;
  lstripPecD1.name = "odd-lstrip-pec-d1";
  lstripPecD1.extent.set(AxisDirection::AxisZ, 1500., 1650.);
  ProtoVolume lstripPecD2;
  lstripPecD2.name = "odd-lstrip-pec-d2";
  lstripPecD2.extent.set(AxisDirection::AxisZ, 1800., 1950.);
  ProtoVolume lstripPecD3;
  lstripPecD3.name = "odd-lstrip-pec-d3";
  lstripPecD3.extent.set(AxisDirection::AxisZ, 2150., 2300.);
  ProtoVolume lstripPecD4;
  lstripPecD4.name = "odd-lstrip-pec-d4";
  lstripPecD4.extent.set(AxisDirection::AxisZ, 2500., 2650.);
  ProtoVolume lstripPecD5;
  lstripPecD5.name = "odd-lstrip-pec-d5";
  lstripPecD5.extent.set(AxisDirection::AxisZ, 2900., 3050.);

  lstripPec.container = ProtoVolume::ContainerStructure{
      {lstripPecD0, lstripPecD1, lstripPecD2, lstripPecD3, lstripPecD4,
       lstripPecD5},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};
  for (auto& cv : lstripPec.container.value().constituentVolumes) {
    cv.internal = ProtoVolume::InternalStructure{Surface::SurfaceType::Disc};
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  lstripContainer.container = ProtoVolume::ContainerStructure{
      {lstripNec, lstripBarrel, lstripPec},
      {BinningData(open, AxisDirection::AxisZ,
                   {-3100., -1200., 1200., 3100.})}};

  // The overall container
  ProtoVolume detectorContainer;
  detectorContainer.name = "odd-light-world";
  detectorContainer.extent.set(AxisDirection::AxisR, 0., 1100.);
  detectorContainer.extent.set(AxisDirection::AxisZ, -3100., 3100.);
  detectorContainer.container = ProtoVolume::ContainerStructure{
      {beamPipeContainer, pixelContainer, pstContainer, sstripContainer,
       lstripContainer},
      {BinningData(open, AxisDirection::AxisR,
                   {0., 25., 200., 210., 720., 1100.})}};

  // ----------------------------------------------------------
  ProtoDetector detector;
  detector.name = "odd-light";
  detector.worldVolume = detectorContainer;

  // Transform into json
  nlohmann::json jdet;
  jdet["detector"] = detector;

  std::ofstream out;
  out.open("odd-proto-detector.json");
  out << jdet.dump(4);
  out.close();

  ProtoDetector detectorIn = jdet["detector"];

  // Let's compare
  BOOST_CHECK_EQUAL(detector.name, detectorIn.name);

  const auto& world = detector.worldVolume;
  const auto& worldIn = detectorIn.worldVolume;

  BOOST_CHECK(isEqual(world, worldIn, 0.1));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
