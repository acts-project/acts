// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE CuboidVolumeBuilderTest
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include <vector>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_CASE(CuboidVolumeBuilderTest)
  {

    // Construct builder
    CuboidVolumeBuilder cvb;

    // Create configurations for surfaces
    std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
    for (unsigned int i = 1; i < 5; i++) {
      // Position of the surfaces
      CuboidVolumeBuilder::SurfaceConfig cfg;
      cfg.position = {i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;

      cfg.detElementConstructor
          = [](std::shared_ptr<const Transform3D>     trans,
               std::shared_ptr<const RectangleBounds> bounds,
               double                                 thickness) {
              return new DetectorElementStub(trans, bounds, thickness);
            };
      surfaceConfig.push_back(cfg);
    }

    // Test that there are actually 4 surface configurations
    BOOST_CHECK_EQUAL(surfaceConfig.size(), 4);

    // Test that 4 surfaces can be built
    for (const auto& cfg : surfaceConfig) {
      std::shared_ptr<const PlaneSurface> pSur = cvb.buildSurface(cfg);
      BOOST_CHECK_NE(pSur, nullptr);
      CHECK_CLOSE_ABS(pSur->center(), cfg.position, 1e-9);
      BOOST_CHECK_NE(pSur->associatedMaterial(), nullptr);
      BOOST_CHECK_NE(pSur->associatedDetectorElement(), nullptr);
    }

    ////////////////////////////////////////////////////////////////////
    // Build layer configurations
    std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig;
    for (auto& sCfg : surfaceConfig) {
      CuboidVolumeBuilder::LayerConfig cfg;
      cfg.surfaceCfg = sCfg;
      layerConfig.push_back(cfg);
    }

    // Test that there are actually 4 layer configurations
    BOOST_CHECK_EQUAL(layerConfig.size(), 4);

    // Test that 4 layers with surfaces can be built
    for (auto& cfg : layerConfig) {
      LayerPtr layer = cvb.buildLayer(cfg);
      BOOST_CHECK_NE(layer, nullptr);
      BOOST_CHECK_NE(cfg.surface, nullptr);
      BOOST_CHECK_EQUAL(layer->surfaceArray()->surfaces().size(), 1);
      BOOST_CHECK_EQUAL(layer->layerType(), LayerType::active);
    }

    for (auto& cfg : layerConfig) {
      cfg.surface = nullptr;
    }

    // Build volume configuration
    CuboidVolumeBuilder::VolumeConfig volumeConfig;
    volumeConfig.position = {2.5 * units::_m, 0., 0.};
    volumeConfig.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig.layerCfg = layerConfig;
    volumeConfig.name     = "Test volume";
    volumeConfig.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));

    // Test the building
    std::shared_ptr<TrackingVolume> trVol = cvb.buildVolume(volumeConfig);
    BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4);
    BOOST_CHECK_EQUAL(trVol->confinedLayers()->arrayObjects().size(),
                      volumeConfig.layers.size() * 2
                          + 1);  // #layers = navigation + material layers
    BOOST_CHECK_EQUAL(trVol->volumeName(), volumeConfig.name);
    BOOST_CHECK_NE(trVol->material(), nullptr);

    // Test the building
    volumeConfig.layers.clear();
    trVol = cvb.buildVolume(volumeConfig);
    BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4);
    BOOST_CHECK_EQUAL(trVol->confinedLayers()->arrayObjects().size(),
                      volumeConfig.layers.size() * 2
                          + 1);  // #layers = navigation + material layers
    BOOST_CHECK_EQUAL(trVol->volumeName(), volumeConfig.name);

    volumeConfig.layers.clear();
    for (auto& lay : volumeConfig.layerCfg) {
      lay.surface = nullptr;
      lay.active  = true;
    }
    trVol = cvb.buildVolume(volumeConfig);
    BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4);
    for (auto& lay : volumeConfig.layers) {
      BOOST_CHECK_EQUAL(lay->layerType(), LayerType::active);
    }

    volumeConfig.layers.clear();
    for (auto& lay : volumeConfig.layerCfg) {
      lay.active = true;
    }
    trVol = cvb.buildVolume(volumeConfig);
    BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4);
    for (auto& lay : volumeConfig.layers) {
      BOOST_CHECK_EQUAL(lay->layerType(), LayerType::active);
    }

    ////////////////////////////////////////////////////////////////////
    // Build TrackingGeometry configuration

    // Build second volume
    std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig2;
    for (int i = 1; i < 5; i++) {
      // Position of the surfaces
      CuboidVolumeBuilder::SurfaceConfig cfg;
      cfg.position = {-i * units::_m, 0., 0.};

      // Rotation of the surfaces
      double   rotationAngle = M_PI * 0.5;
      Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
      Vector3D yPos(0., 1., 0.);
      Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
      cfg.rotation.col(0) = xPos;
      cfg.rotation.col(1) = yPos;
      cfg.rotation.col(2) = zPos;

      // Boundaries of the surfaces
      cfg.rBounds = std::make_shared<const RectangleBounds>(
          RectangleBounds(0.5 * units::_m, 0.5 * units::_m));

      // Material of the surfaces
      MaterialProperties matProp(
          352.8, 407., 9.012, 4., 1.848e-3, 0.5 * units::_mm);
      cfg.surMat = std::shared_ptr<const SurfaceMaterial>(
          new HomogeneousSurfaceMaterial(matProp));

      // Thickness of the detector element
      cfg.thickness = 1. * units::_um;
      surfaceConfig2.push_back(cfg);
    }

    std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig2;
    for (auto& sCfg : surfaceConfig2) {
      CuboidVolumeBuilder::LayerConfig cfg;
      cfg.surfaceCfg = sCfg;
      layerConfig2.push_back(cfg);
    }
    CuboidVolumeBuilder::VolumeConfig volumeConfig2;
    volumeConfig2.position = {-2.5 * units::_m, 0., 0.};
    volumeConfig2.length   = {5. * units::_m, 1. * units::_m, 1. * units::_m};
    volumeConfig2.layerCfg = layerConfig2;
    volumeConfig2.name     = "Test volume2";

    CuboidVolumeBuilder::Config config;
    config.position  = {0., 0., 0.};
    config.length    = {10. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {volumeConfig2, volumeConfig};

    cvb.setConfig(config);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder tgb(tgbCfg);

    std::unique_ptr<const TrackingGeometry> detector = tgb.trackingGeometry();
    BOOST_CHECK_EQUAL(
        detector->lowestTrackingVolume({1., 0., 0.})->volumeName(),
        volumeConfig.name);
    BOOST_CHECK_EQUAL(
        detector->lowestTrackingVolume({-1., 0., 0.})->volumeName(),
        volumeConfig2.name);
  }
}  // namespace Test
}  // namespace Acts
