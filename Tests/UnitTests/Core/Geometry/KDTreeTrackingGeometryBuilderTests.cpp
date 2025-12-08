// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(KDTreeTrackingGeometryBuilder_simple) {
  GeometryContext tContext;
  CylindricalTrackingGeometry ctGeometry(tContext);
  CylindricalTrackingGeometry::DetectorStore detectorStore;

  // The collected surfaces
  std::vector<std::shared_ptr<Surface>> layerSurfacePtrs;

  // Add a beam pipe
  auto hTransform = Transform3::Identity();
  layerSurfacePtrs.push_back(
      Surface::makeShared<CylinderSurface>(hTransform, 15., 800.));

  // Pixel Surfaces
  std::vector<double> pLayerRadii = {32., 72., 116., 172.};
  std::vector<std::pair<int, int>> pLayerBinning = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  std::vector<double> pModuleTiltPhi = {0.145, 0.145, 0.145, 0.145};
  std::vector<double> pModuleHalfX = {8.4, 8.4, 8.4, 8.4};
  std::vector<double> pModuleHalfY = {36., 36., 36., 36.};
  std::vector<double> pModuleThickness = {0.15, 0.15, 0.15, 0.15};

  // Fill surfaces from cylinder layers
  for (std::size_t ilp = 0; ilp < pLayerRadii.size(); ++ilp) {
    std::vector<Surface*> layerSurfaces = ctGeometry.surfacesCylinder(
        detectorStore, pModuleHalfX[ilp], pModuleHalfY[ilp],
        pModuleThickness[ilp], pModuleTiltPhi[ilp], pLayerRadii[ilp], 2_mm,
        5_mm, pLayerBinning[ilp]);

    // Make a shared version out of it
    for (auto& sf : layerSurfaces) {
      layerSurfacePtrs.push_back(sf->getSharedPtr());
    }
  }

  // Fill surfaces for disc layers
  std::vector<double> discZ = {-700., -600., 600., 700.};
  std::vector<double> discRadii = {60., 60., 60., 60.};
  std::vector<int> discModules = {22, 22, 22, 22};

  std::vector<double> dModuleHalfXMinY = {6.4, 6.4, 6.4, 6.4};
  std::vector<double> dModuleHalfXMaxY = {12.4, 12.4, 12.4, 12.4};
  std::vector<double> dModuleHalfY = {36., 36., 36., 36.};
  std::vector<double> dModuleTilt = {0.075, 0.075, 0.075, 0.075};
  std::vector<double> dModuleThickness = {0.15, 0.15, 0.15, 0.15};

  for (std::size_t ilp = 0; ilp < discZ.size(); ++ilp) {
    std::vector<Surface*> layerSurfaces = ctGeometry.surfacesRing(
        detectorStore, dModuleHalfXMinY[ilp], dModuleHalfXMaxY[ilp],
        dModuleHalfY[ilp], dModuleThickness[ilp], dModuleTilt[ilp],
        discRadii[ilp], discZ[ilp], 2., discModules[ilp]);
    for (auto& sf : layerSurfaces) {
      layerSurfacePtrs.push_back(sf->getSharedPtr());
    }
  }

  // Make a proto detectpr description
  ProtoVolume beamPipeContainer;
  beamPipeContainer.name = "odd-beam-pipe";
  beamPipeContainer.extent.set(AxisDirection::AxisR, 0., 17);
  ProtoVolume beamPipe;
  beamPipe.name = "odd-beam-pipe-l";
  beamPipe.extent.set(AxisDirection::AxisR, 2., 16.);
  beamPipe.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  beamPipeContainer.container = ProtoVolume::ContainerStructure{
      {beamPipe}, {BinningData(open, AxisDirection::AxisR, {0., 1.})}, true};

  // Pixel section
  ProtoVolume pixelContainer;
  pixelContainer.name = "odd-pixel";
  pixelContainer.extent.set(AxisDirection::AxisR, 18., 200);

  ProtoVolume pixelNec;
  pixelNec.name = "odd-pixel-nec";
  pixelNec.extent.set(AxisDirection::AxisZ, -1000., -580);

  ProtoVolume pixNecD1;
  pixNecD1.name = "odd-pixel-nec-d1";
  pixNecD1.extent.set(AxisDirection::AxisZ, -720., -680);
  ProtoVolume pixNecD0;
  pixNecD0.name = "odd-pixel-nec-d0";
  pixNecD0.extent.set(AxisDirection::AxisZ, -620., -580);
  pixelNec.container = ProtoVolume::ContainerStructure{
      {pixNecD1, pixNecD0},
      {BinningData(open, AxisDirection::AxisZ, {0., 1.})},
      true};
  for (auto& cv : pixelNec.container.value().constituentVolumes) {
    cv.internal = ProtoVolume::InternalStructure{Surface::SurfaceType::Disc};
  }

  ProtoVolume pixelBarrel;
  pixelBarrel.name = "odd-pixel-barrel";
  pixelBarrel.extent.set(AxisDirection::AxisZ, -580., 580);

  ProtoVolume pixBarrelL0;
  pixBarrelL0.name = "odd-pixel-barrel-l0";
  pixBarrelL0.extent.set(AxisDirection::AxisR, 28., 48.);
  pixBarrelL0.extent.set(AxisDirection::AxisZ, -580., 580);
  pixBarrelL0.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  ProtoVolume pixBarrelL1;
  pixBarrelL1.name = "odd-pixel-barrel-l1";
  pixBarrelL1.extent.set(AxisDirection::AxisR, 62., 76);
  pixBarrelL1.extent.set(AxisDirection::AxisZ, -580., 580);
  pixBarrelL1.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  ProtoVolume pixBarrelL2;
  pixBarrelL2.name = "odd-pixel-barrel-l2";
  pixBarrelL2.extent.set(AxisDirection::AxisR, 100., 120.);
  pixBarrelL2.extent.set(AxisDirection::AxisZ, -580., 580);
  pixBarrelL2.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};
  ProtoVolume pixBarrelL3;
  pixBarrelL3.name = "odd-pixel-barrel-l3";
  pixBarrelL3.extent.set(AxisDirection::AxisR, 160., 180.);
  pixBarrelL3.extent.set(AxisDirection::AxisZ, -580., 580);
  pixBarrelL3.internal =
      ProtoVolume::InternalStructure{Surface::SurfaceType::Cylinder};

  pixelBarrel.container = ProtoVolume::ContainerStructure{
      {pixBarrelL0, pixBarrelL1, pixBarrelL2, pixBarrelL3},
      {BinningData(open, AxisDirection::AxisR, {0., 1})},
      true};

  ProtoVolume pixelPec;
  pixelPec.name = "odd-pixel-pec";
  pixelPec.extent.set(AxisDirection::AxisZ, 580., 1000.);

  ProtoVolume pixPecD0;
  pixPecD0.name = "odd-pixel-pec-d0";
  pixPecD0.extent.set(AxisDirection::AxisZ, 580., 620);
  ProtoVolume pixPecD1;
  pixPecD1.name = "odd-pixel-pec-d1";
  pixPecD1.extent.set(AxisDirection::AxisZ, 680., 720);

  pixelPec.container = ProtoVolume::ContainerStructure{
      {pixPecD0, pixPecD1},
      {BinningData(open, AxisDirection::AxisZ, {0., 1})},
      true};
  for (auto& cv : pixelPec.container.value().constituentVolumes) {
    cv.internal = ProtoVolume::InternalStructure{Surface::SurfaceType::Disc};
  }

  pixelContainer.container = ProtoVolume::ContainerStructure{
      {pixelNec, pixelBarrel, pixelPec},
      {BinningData(open, AxisDirection::AxisZ, {-1000., -580., 580., 1000.})}};

  ProtoVolume detectorContainer;
  detectorContainer.name = "odd-detector";
  detectorContainer.extent.set(AxisDirection::AxisR, 0., 200);
  detectorContainer.container = ProtoVolume::ContainerStructure{
      {beamPipeContainer, pixelContainer},
      {BinningData(open, AxisDirection::AxisR, {0., 17.5, 200.})}};

  ProtoDetector detector;
  detector.name = "odd";
  detector.worldVolume = detectorContainer;

  auto logLevel = Logging::VERBOSE;

  // Surface array creator
  auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
      SurfaceArrayCreator::Config(),
      getDefaultLogger("SurfaceArrayCreator", logLevel));
  // Layer Creator
  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<LayerCreator>(
      lcConfig, getDefaultLogger("LayerCreator", logLevel));
  // Layer array creator
  LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
      lacConfig, getDefaultLogger("LayerArrayCreator", logLevel));
  // Tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      tvacConfig, getDefaultLogger("TrackingVolumeArrayCreator", logLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", logLevel));

  // The KDT tracking geometry builder
  KDTreeTrackingGeometryBuilder::Config kdtgConfig;
  kdtgConfig.layerCreator = layerCreator;
  kdtgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  // Reserve the right amount of surfaces
  kdtgConfig.surfaces = layerSurfacePtrs;

  // Assign the proto detector
  kdtgConfig.protoDetector = detector;

  // Make the builder
  auto kdtTrackingGeometryBuilder = KDTreeTrackingGeometryBuilder(
      kdtgConfig, getDefaultLogger("KDTreeTrackingGeometryBuilder", logLevel));

  auto trackingGeometry = kdtTrackingGeometryBuilder.trackingGeometry(tContext);
  BOOST_CHECK(trackingGeometry != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
