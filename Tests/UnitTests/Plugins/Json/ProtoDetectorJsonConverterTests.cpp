// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/ProtoDetectorJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <fstream>
#include <iostream>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(ProtoDetectorJsonConverter)

BOOST_AUTO_TEST_CASE(ProtoDetectorRoundTrip) {
  Acts::ExtentEnvelope cylinderLayerEnvelope = zeroEnvelopes;
  cylinderLayerEnvelope[Acts::binR] = {1., 1.};
  cylinderLayerEnvelope[Acts::binZ] = {2., 2.};

  Acts::ExtentEnvelope discLayerEnvelope = zeroEnvelopes;
  discLayerEnvelope[Acts::binR] = {1., 1.};
  discLayerEnvelope[Acts::binZ] = {1., 1.};

  // Beam Pipe container
  Acts::ProtoVolume beamPipeContainer;
  beamPipeContainer.name = "odd-beam-pipe";
  beamPipeContainer.extent.set(Acts::binR, 0., 25);
  Acts::ProtoVolume beamPipe;
  beamPipe.name = "odd-beam-pipe-l";
  beamPipe.extent.set(Acts::binR, 2., 24.);
  beamPipe.layerType = Acts::Surface::SurfaceType::Cylinder;
  beamPipeContainer.constituentVolumes = {beamPipe};
  beamPipeContainer.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1.})};
  beamPipeContainer.layerContainer = true;

  // Pixel section
  Acts::ProtoVolume pixelContainer;
  pixelContainer.name = "odd-pixel";
  pixelContainer.extent.set(Acts::binR, 25., 200);

  Acts::ProtoVolume pixelNec;
  pixelNec.name = "odd-pixel-nec";
  pixelNec.extent.set(Acts::binZ, -3100., -580);

  Acts::ProtoVolume pixNecD6;
  pixNecD6.name = "odd-pixel-nec-d6";
  pixNecD6.extent.set(Acts::binZ, -1540., -1500);
  Acts::ProtoVolume pixNecD5;
  pixNecD5.name = "odd-pixel-nec-d5";
  pixNecD5.extent.set(Acts::binZ, -1340., -1300);
  Acts::ProtoVolume pixNecD4;
  pixNecD4.name = "odd-pixel-nec-d4";
  pixNecD4.extent.set(Acts::binZ, -1140., -1100);
  Acts::ProtoVolume pixNecD3;
  pixNecD3.name = "odd-pixel-nec-d3";
  pixNecD3.extent.set(Acts::binZ, -1000., -960.);
  Acts::ProtoVolume pixNecD2;
  pixNecD2.name = "odd-pixel-nec-d2";
  pixNecD2.extent.set(Acts::binZ, -860., -820);
  Acts::ProtoVolume pixNecD1;
  pixNecD1.name = "odd-pixel-nec-d1";
  pixNecD1.extent.set(Acts::binZ, -740., -700);
  Acts::ProtoVolume pixNecD0;
  pixNecD0.name = "odd-pixel-nec-d0";
  pixNecD0.extent.set(Acts::binZ, -640., -600);
  pixelNec.constituentVolumes = {pixNecD6, pixNecD5, pixNecD4, pixNecD3,
                                 pixNecD2, pixNecD1, pixNecD0};

  Acts::BinningData pixEcBinningR =
      Acts::BinningData(Acts::open, Acts::binR, 2., 0., 1.);
  Acts::BinningData pixEcBinningPhi =
      Acts::BinningData(Acts::closed, Acts::binPhi, 30., -M_PI, M_PI);

  pixelNec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1.})};
  for (auto& cv : pixelNec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.layerSurfaceBinning = {pixEcBinningR, pixEcBinningPhi};
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  pixelNec.layerContainer = true;

  Acts::ProtoVolume pixelBarrel;
  pixelBarrel.name = "odd-pixel-barrel";
  pixelBarrel.extent.set(Acts::binZ, -580., 580);

  Acts::ProtoVolume pixBarrelL0;
  pixBarrelL0.name = "odd-pixel-barrel-l0";
  pixBarrelL0.extent.set(Acts::binR, 28., 48.);
  pixBarrelL0.extent.set(Acts::binZ, -580., 580);
  Acts::ProtoVolume pixBarrelL1;
  pixBarrelL1.name = "odd-pixel-barrel-l1";
  pixBarrelL1.extent.set(Acts::binR, 62., 76);
  pixBarrelL1.extent.set(Acts::binZ, -580., 580);
  Acts::ProtoVolume pixBarrelL2;
  pixBarrelL2.name = "odd-pixel-barrel-l2";
  pixBarrelL2.extent.set(Acts::binR, 100., 120.);
  pixBarrelL2.extent.set(Acts::binZ, -580., 580);
  Acts::ProtoVolume pixBarrelL3;
  pixBarrelL3.name = "odd-pixel-barrel-l3";
  pixBarrelL3.extent.set(Acts::binR, 160., 180.);
  pixBarrelL3.extent.set(Acts::binZ, -580., 580);

  pixelBarrel.constituentVolumes = {pixBarrelL0, pixBarrelL1, pixBarrelL2,
                                    pixBarrelL3};
  for (auto& cv : pixelBarrel.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Cylinder;
    cv.extent.setEnvelope(cylinderLayerEnvelope);
  }

  pixelBarrel.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1})};
  pixelBarrel.layerContainer = true;

  Acts::ProtoVolume pixelPec;
  pixelPec.name = "odd-pixel-pec";
  pixelPec.extent.set(Acts::binZ, 580., 3100.);

  Acts::ProtoVolume pixPecD0;
  pixPecD0.name = "odd-pixel-pec-d0";
  pixPecD0.extent.set(Acts::binZ, 600., 640);
  Acts::ProtoVolume pixPecD1;
  pixPecD1.name = "odd-pixel-pec-d1";
  pixPecD1.extent.set(Acts::binZ, 700., 740);
  Acts::ProtoVolume pixPecD2;
  pixPecD2.name = "odd-pixel-pec-d2";
  pixPecD2.extent.set(Acts::binZ, 820., 860.);
  Acts::ProtoVolume pixPecD3;
  pixPecD3.name = "odd-pixel-pec-d3";
  pixPecD3.extent.set(Acts::binZ, 960., 1000.);
  Acts::ProtoVolume pixPecD4;
  pixPecD4.name = "odd-pixel-pec-d4";
  pixPecD4.extent.set(Acts::binZ, 1100., 1140);
  Acts::ProtoVolume pixPecD5;
  pixPecD5.name = "odd-pixel-pec-d5";
  pixPecD5.extent.set(Acts::binZ, 1300., 1340.);
  Acts::ProtoVolume pixPecD6;
  pixPecD6.name = "odd-pixel-pec-d6";
  pixPecD6.extent.set(Acts::binZ, 1500., 1540.);

  pixelPec.constituentVolumes = {pixPecD0, pixPecD1, pixPecD2, pixPecD3,
                                 pixPecD4, pixPecD5, pixPecD6};
  pixelPec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};
  for (auto& cv : pixelPec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.layerSurfaceBinning = {pixEcBinningR, pixEcBinningPhi};
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  pixelPec.layerContainer = true;

  pixelContainer.constituentVolumes = {pixelNec, pixelBarrel, pixelPec};
  pixelContainer.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {-3100., -580., 580., 3100.})};

  // Short Strip section
  Acts::ProtoVolume pstContainer;
  pstContainer.name = "odd-pst";
  pstContainer.extent.set(Acts::binR, 200., 210.);
  Acts::ProtoVolume pst;
  pst.name = "odd-pst-l";
  pst.extent.set(Acts::binR, 201., 209.);
  pst.layerType = Acts::Surface::SurfaceType::Cylinder;
  pstContainer.constituentVolumes = {pst};
  pstContainer.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1.})};
  pstContainer.layerContainer = true;

  // Short Strip section
  Acts::ProtoVolume sstripContainer;
  sstripContainer.name = "odd-sstrip";
  sstripContainer.extent.set(Acts::binR, 210., 720);

  Acts::BinningData sstripEcBinningR =
      Acts::BinningData(Acts::open, Acts::binR, 3., 0., 1.);
  Acts::BinningData sstripEcBinningPhi =
      Acts::BinningData(Acts::closed, Acts::binPhi, 42., -M_PI, M_PI);

  Acts::ProtoVolume sstripNec;
  sstripNec.name = "odd-sstrip-nec";
  sstripNec.extent.set(Acts::binZ, -3100., -1200);
  Acts::ProtoVolume sstripNecD5;
  sstripNecD5.name = "odd-sstrip-nec-d5";
  sstripNecD5.extent.set(Acts::binZ, -3000, -2900.);
  Acts::ProtoVolume sstripNecD4;
  sstripNecD4.name = "odd-sstrip-nec-d4";
  sstripNecD4.extent.set(Acts::binZ, -2600., -2500.);
  Acts::ProtoVolume sstripNecD3;
  sstripNecD3.name = "odd-sstrip-nec-d3";
  sstripNecD3.extent.set(Acts::binZ, -2250, -2150.);
  Acts::ProtoVolume sstripNecD2;
  sstripNecD2.name = "odd-sstrip-nec-d2";
  sstripNecD2.extent.set(Acts::binZ, -1900, -1800.);
  Acts::ProtoVolume sstripNecD1;
  sstripNecD1.name = "odd-sstrip-nec-d1";
  sstripNecD1.extent.set(Acts::binZ, -1600., -1500.);
  Acts::ProtoVolume sstripNecD0;
  sstripNecD0.name = "odd-sstrip-nec-d0";
  sstripNecD0.extent.set(Acts::binZ, -1350., -1250.);

  sstripNec.constituentVolumes = {sstripNecD5, sstripNecD4, sstripNecD3,
                                  sstripNecD2, sstripNecD1, sstripNecD0};
  sstripNec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};
  for (auto& cv : sstripNec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.layerSurfaceBinning = {sstripEcBinningR, sstripEcBinningPhi};
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  sstripNec.layerContainer = true;

  Acts::ProtoVolume sstripBarrel;
  sstripBarrel.name = "odd-sstrip-barrel";
  sstripBarrel.extent.set(Acts::binZ, -1200., 1200);

  Acts::ProtoVolume sstripBarrelL0;
  sstripBarrelL0.name = "odd-sstrip-barrel-l0";
  sstripBarrelL0.extent.set(Acts::binR, 240., 280.);
  Acts::ProtoVolume sstripBarrelL1;
  sstripBarrelL1.name = "odd-sstrip-barrel-l1";
  sstripBarrelL1.extent.set(Acts::binR, 340., 380.);
  Acts::ProtoVolume sstripBarrelL2;
  sstripBarrelL2.name = "odd-sstrip-barrel-l2";
  sstripBarrelL2.extent.set(Acts::binR, 480., 520.);
  Acts::ProtoVolume sstripBarrelL3;
  sstripBarrelL3.name = "odd-sstrip-barrel-l3";
  sstripBarrelL3.extent.set(Acts::binR, 640., 680.);

  sstripBarrel.constituentVolumes = {sstripBarrelL0, sstripBarrelL1,
                                     sstripBarrelL2, sstripBarrelL3};
  for (auto& cv : sstripBarrel.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Cylinder;
    cv.extent.setEnvelope(cylinderLayerEnvelope);
  }

  sstripBarrel.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1})};
  sstripBarrel.layerContainer = true;

  Acts::ProtoVolume sstripPec;
  sstripPec.name = "odd-sstrip-pec";
  sstripPec.extent.set(Acts::binZ, 1200., 3100);

  Acts::ProtoVolume sstripPecD0;
  sstripPecD0.name = "odd-sstrip-pec-d0";
  sstripPecD0.extent.set(Acts::binZ, 1250., 1350);
  Acts::ProtoVolume sstripPecD1;
  sstripPecD1.name = "odd-sstrip-pec-d1";
  sstripPecD1.extent.set(Acts::binZ, 1500., 1600.);
  Acts::ProtoVolume sstripPecD2;
  sstripPecD2.name = "odd-sstrip-pec-d2";
  sstripPecD2.extent.set(Acts::binZ, 1800., 1900.);
  Acts::ProtoVolume sstripPecD3;
  sstripPecD3.name = "odd-sstrip-pec-d3";
  sstripPecD3.extent.set(Acts::binZ, 2150., 2250.);
  Acts::ProtoVolume sstripPecD4;
  sstripPecD4.name = "odd-sstrip-pec-d4";
  sstripPecD4.extent.set(Acts::binZ, 2500., 2600.);
  Acts::ProtoVolume sstripPecD5;
  sstripPecD5.name = "odd-sstrip-pec-d5";
  sstripPecD5.extent.set(Acts::binZ, 2900., 3000.);

  sstripPec.constituentVolumes = {sstripPecD0, sstripPecD1, sstripPecD2,
                                  sstripPecD3, sstripPecD4, sstripPecD5};
  sstripPec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};
  for (auto& cv : sstripPec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.layerSurfaceBinning = {sstripEcBinningR, sstripEcBinningPhi};
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  sstripPec.layerContainer = true;

  sstripContainer.constituentBinning = {Acts::BinningData(
      Acts::open, Acts::binZ, {-3100., -1200., 1200., 3100.})};
  sstripContainer.constituentVolumes = {sstripNec, sstripBarrel, sstripPec};

  // Long Strip section
  Acts::ProtoVolume lstripContainer;
  lstripContainer.name = "odd-lstrip";
  lstripContainer.extent.set(Acts::binR, 720, 1100.);

  Acts::ProtoVolume lstripNec;
  lstripNec.name = "odd-lstrip-nec";
  lstripNec.extent.set(Acts::binZ, -3100., -1200);
  Acts::ProtoVolume lstripNecD5;
  lstripNecD5.name = "odd-lstrip-nec-d5";
  lstripNecD5.extent.set(Acts::binZ, -3050, -2900.);
  Acts::ProtoVolume lstripNecD4;
  lstripNecD4.name = "odd-lstrip-nec-d4";
  lstripNecD4.extent.set(Acts::binZ, -2650., -2500.);
  Acts::ProtoVolume lstripNecD3;
  lstripNecD3.name = "odd-lstrip-nec-d3";
  lstripNecD3.extent.set(Acts::binZ, -2300, -2150.);
  Acts::ProtoVolume lstripNecD2;
  lstripNecD2.name = "odd-lstrip-nec-d2";
  lstripNecD2.extent.set(Acts::binZ, -1950, -1800.);
  Acts::ProtoVolume lstripNecD1;
  lstripNecD1.name = "odd-lstrip-nec-d1";
  lstripNecD1.extent.set(Acts::binZ, -1650., -1500.);
  Acts::ProtoVolume lstripNecD0;
  lstripNecD0.name = "odd-lstrip-nec-d0";
  lstripNecD0.extent.set(Acts::binZ, -1400., -1250.);

  lstripNec.constituentVolumes = {lstripNecD5, lstripNecD4, lstripNecD3,
                                  lstripNecD2, lstripNecD1, lstripNecD0};
  lstripNec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};
  for (auto& cv : lstripNec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  lstripNec.layerContainer = true;

  Acts::ProtoVolume lstripBarrel;
  lstripBarrel.name = "odd-lstrip-barrel";
  lstripBarrel.extent.set(Acts::binZ, -1200., 1200);

  Acts::ProtoVolume lstripBarrelL0;
  lstripBarrelL0.name = "odd-lstrip-barrel-l0";
  lstripBarrelL0.extent.set(Acts::binR, 800., 840.);
  Acts::ProtoVolume lstripBarrelL1;
  lstripBarrelL1.name = "odd-lstrip-barrel-l1";
  lstripBarrelL1.extent.set(Acts::binR, 1000., 1050.);

  lstripBarrel.constituentVolumes = {lstripBarrelL0, lstripBarrelL1};
  for (auto& cv : lstripBarrel.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Cylinder;
    cv.extent.setEnvelope(cylinderLayerEnvelope);
  }

  lstripBarrel.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binR, {0., 1})};
  lstripBarrel.layerContainer = true;

  Acts::ProtoVolume lstripPec;
  lstripPec.name = "odd-lstrip-pec";
  lstripPec.extent.set(Acts::binZ, 1200., 3100);

  Acts::ProtoVolume lstripPecD0;
  lstripPecD0.name = "odd-lstrip-pec-d0";
  lstripPecD0.extent.set(Acts::binZ, 1250., 1400);
  Acts::ProtoVolume lstripPecD1;
  lstripPecD1.name = "odd-lstrip-pec-d1";
  lstripPecD1.extent.set(Acts::binZ, 1500., 1650.);
  Acts::ProtoVolume lstripPecD2;
  lstripPecD2.name = "odd-lstrip-pec-d2";
  lstripPecD2.extent.set(Acts::binZ, 1800., 1950.);
  Acts::ProtoVolume lstripPecD3;
  lstripPecD3.name = "odd-lstrip-pec-d3";
  lstripPecD3.extent.set(Acts::binZ, 2150., 2300.);
  Acts::ProtoVolume lstripPecD4;
  lstripPecD4.name = "odd-lstrip-pec-d4";
  lstripPecD4.extent.set(Acts::binZ, 2500., 2650.);
  Acts::ProtoVolume lstripPecD5;
  lstripPecD5.name = "odd-lstrip-pec-d5";
  lstripPecD5.extent.set(Acts::binZ, 2900., 3050.);

  lstripPec.constituentVolumes = {lstripPecD0, lstripPecD1, lstripPecD2,
                                  lstripPecD3, lstripPecD4, lstripPecD5};
  lstripPec.constituentBinning = {
      Acts::BinningData(Acts::open, Acts::binZ, {0., 1})};
  for (auto& cv : lstripPec.constituentVolumes) {
    cv.layerType = Acts::Surface::SurfaceType::Disc;
    cv.extent.setEnvelope(discLayerEnvelope);
  }
  lstripPec.layerContainer = true;

  lstripContainer.constituentVolumes = {lstripNec, lstripBarrel, lstripPec};
  lstripContainer.constituentBinning = {Acts::BinningData(
      Acts::open, Acts::binZ, {-3100., -1200., 1200., 3100.})};

  // The overall container
  Acts::ProtoVolume detectorContainer;
  detectorContainer.name = "odd-light-world";
  detectorContainer.extent.set(Acts::binR, 0., 1100.);
  detectorContainer.extent.set(Acts::binZ, -3100., 3100.);
  detectorContainer.constituentVolumes = {beamPipeContainer, pixelContainer,
                                          pstContainer, sstripContainer,
                                          lstripContainer};
  detectorContainer.constituentBinning = {Acts::BinningData(
      Acts::open, Acts::binR, {0., 25., 200., 210., 720., 1100.})};

  // ----------------------------------------------------------
  Acts::ProtoDetector detector;
  detector.name = "odd-light";
  detector.worldVolume = detectorContainer;

  // Transform into json
  nlohmann::json jdet;
  jdet["detector"] = detector;

  std::ofstream out;
  out.open("odd-proto-detector.json");
  out << jdet.dump(4);
  out.close();

  Acts::ProtoDetector detectorIn = jdet["detector"];
}

BOOST_AUTO_TEST_SUITE_END()
