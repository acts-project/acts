// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE CuboidVolumeBuilderTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <vector>

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_CASE(CuboidVolumeBuilderTest) {
  // Construct builder
  CuboidVolumeBuilder cvb;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Create configurations for surfaces
  std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
  for (unsigned int i = 1; i < 5; i++) {
    // Position of the surfaces
    CuboidVolumeBuilder::SurfaceConfig cfg;
    cfg.position = {i * UnitConstants::m, 0., 0.};

    // Rotation of the surfaces
    double rotationAngle = M_PI * 0.5;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    cfg.rotation.col(0) = xPos;
    cfg.rotation.col(1) = yPos;
    cfg.rotation.col(2) = zPos;

    // Boundaries of the surfaces
    cfg.rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));

    // Material of the surfaces
    MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5_mm);
    cfg.surMat = std::shared_ptr<const ISurfaceMaterial>(
        new HomogeneousSurfaceMaterial(matProp));

    // Thickness of the detector element
    cfg.thickness = 1_um;

    cfg.detElementConstructor =
        [](std::shared_ptr<const Transform3D> trans,
           std::shared_ptr<const RectangleBounds> bounds, double thickness) {
          return new DetectorElementStub(trans, bounds, thickness);
        };
    surfaceConfig.push_back(cfg);
  }

  // Test that there are actually 4 surface configurations
  BOOST_CHECK_EQUAL(surfaceConfig.size(), 4u);

  // Test that 4 surfaces can be built
  for (const auto& cfg : surfaceConfig) {
    std::shared_ptr<const PlaneSurface> pSur = cvb.buildSurface(tgContext, cfg);
    BOOST_CHECK_NE(pSur, nullptr);
    CHECK_CLOSE_ABS(pSur->center(tgContext), cfg.position, 1e-9);
    BOOST_CHECK_NE(pSur->surfaceMaterial(), nullptr);
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
  BOOST_CHECK_EQUAL(layerConfig.size(), 4u);

  // Test that 4 layers with surfaces can be built
  for (auto& cfg : layerConfig) {
    LayerPtr layer = cvb.buildLayer(tgContext, cfg);
    BOOST_CHECK_NE(layer, nullptr);
    BOOST_CHECK_NE(cfg.surface, nullptr);
    BOOST_CHECK_EQUAL(layer->surfaceArray()->surfaces().size(), 1u);
    BOOST_CHECK_EQUAL(layer->layerType(), LayerType::active);
  }

  for (auto& cfg : layerConfig) {
    cfg.surface = nullptr;
  }

  // Build volume configuration
  CuboidVolumeBuilder::VolumeConfig volumeConfig;
  volumeConfig.position = {2.5_m, 0., 0.};
  volumeConfig.length = {5_m, 1_m, 1_m};
  volumeConfig.layerCfg = layerConfig;
  volumeConfig.name = "Test volume";
  volumeConfig.volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(
          Material(352.8, 407., 9.012, 4., 1.848e-3));

  // Test the building
  std::shared_ptr<TrackingVolume> trVol =
      cvb.buildVolume(tgContext, volumeConfig);
  BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4u);
  BOOST_CHECK_EQUAL(trVol->confinedLayers()->arrayObjects().size(),
                    volumeConfig.layers.size() * 2 +
                        1u);  // #layers = navigation + material layers
  BOOST_CHECK_EQUAL(trVol->volumeName(), volumeConfig.name);
  BOOST_CHECK_NE(trVol->volumeMaterial(), nullptr);

  // Test the building
  volumeConfig.layers.clear();
  trVol = cvb.buildVolume(tgContext, volumeConfig);
  BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4u);
  BOOST_CHECK_EQUAL(trVol->confinedLayers()->arrayObjects().size(),
                    volumeConfig.layers.size() * 2 +
                        1u);  // #layers = navigation + material layers
  BOOST_CHECK_EQUAL(trVol->volumeName(), volumeConfig.name);

  volumeConfig.layers.clear();
  for (auto& lay : volumeConfig.layerCfg) {
    lay.surface = nullptr;
    lay.active = true;
  }
  trVol = cvb.buildVolume(tgContext, volumeConfig);
  BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4u);
  for (auto& lay : volumeConfig.layers) {
    BOOST_CHECK_EQUAL(lay->layerType(), LayerType::active);
  }

  volumeConfig.layers.clear();
  for (auto& lay : volumeConfig.layerCfg) {
    lay.active = true;
  }
  trVol = cvb.buildVolume(tgContext, volumeConfig);
  BOOST_CHECK_EQUAL(volumeConfig.layers.size(), 4u);
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
    cfg.position = {-i * UnitConstants::m, 0., 0.};

    // Rotation of the surfaces
    double rotationAngle = M_PI * 0.5;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    cfg.rotation.col(0) = xPos;
    cfg.rotation.col(1) = yPos;
    cfg.rotation.col(2) = zPos;

    // Boundaries of the surfaces
    cfg.rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));

    // Material of the surfaces
    MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5_mm);
    cfg.surMat = std::shared_ptr<const ISurfaceMaterial>(
        new HomogeneousSurfaceMaterial(matProp));

    // Thickness of the detector element
    cfg.thickness = 1_um;
    surfaceConfig2.push_back(cfg);
  }

  std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig2;
  for (auto& sCfg : surfaceConfig2) {
    CuboidVolumeBuilder::LayerConfig cfg;
    cfg.surfaceCfg = sCfg;
    layerConfig2.push_back(cfg);
  }
  CuboidVolumeBuilder::VolumeConfig volumeConfig2;
  volumeConfig2.position = {-2.5_m, 0., 0.};
  volumeConfig2.length = {5_m, 1_m, 1_m};
  volumeConfig2.layerCfg = layerConfig2;
  volumeConfig2.name = "Test volume2";

  CuboidVolumeBuilder::Config config;
  config.position = {0., 0., 0.};
  config.length = {10_m, 1_m, 1_m};
  config.volumeCfg = {volumeConfig2, volumeConfig};

  cvb.setConfig(config);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  TrackingGeometryBuilder tgb(tgbCfg);

  std::unique_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);
  BOOST_CHECK_EQUAL(
      detector->lowestTrackingVolume(tgContext, Vector3D(1., 0., 0.))
          ->volumeName(),
      volumeConfig.name);
  BOOST_CHECK_EQUAL(
      detector->lowestTrackingVolume(tgContext, Vector3D(-1., 0., 0.))
          ->volumeName(),
      volumeConfig2.name);
}

  BOOST_AUTO_TEST_CASE(BoxGeometryBuilderTest_confinedVolumes)
  {
    // Production factory
    BoxGeometryBuilder bgb;

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg;
    vCfg.position = {1. * units::_m, 0., 0.};
    vCfg.length   = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg.name     = "Test volume";
    // Build and add 2 confined volumes
    BoxGeometryBuilder::VolumeConfig cvCfg1;
    cvCfg1.position = {1.1 * units::_m, 0., 0.};
    cvCfg1.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg1.name     = "Confined volume1";
    cvCfg1.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    BoxGeometryBuilder::VolumeConfig cvCfg2;
    cvCfg2.position = {0.9 * units::_m, 0., 0.};
    cvCfg2.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg2.name     = "Confined volume2";
    vCfg.volumeCfg  = {cvCfg1, cvCfg2};

    // Build detector
    BoxGeometryBuilder::Config config;
    config.position  = {1. * units::_m, 0., 0.};
    config.length    = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {vCfg};
    std::shared_ptr<TrackingGeometry> detector
        = bgb.buildTrackingGeometry(config);

    // Test that the right volume is selected
    BOOST_TEST(
        detector->lowestTrackingVolume({1. * units::_m, 0., 0.})->volumeName()
        == vCfg.name);
    BOOST_TEST(
        detector->lowestTrackingVolume({1.1 * units::_m, 0., 0.})->volumeName()
        == cvCfg1.name);
    BOOST_TEST(
        detector->lowestTrackingVolume({0.9 * units::_m, 0., 0.})->volumeName()
        == cvCfg2.name);

    // Set propagator and navigator
    PropagatorOptions<ActionList<StepVolumeCollector>> propOpts;
    propOpts.maxStepSize = 10. * units::_mm;
    StraightLineStepper sls;
    Navigator           navi(detector);
    navi.resolvePassive   = true;
    navi.resolveMaterial  = true;
    navi.resolveSensitive = true;

    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        nullptr, startParams, startMom, 1.);

    // Launch and collect results
    const auto& result = prop.propagate(sbtp, propOpts);
    const StepVolumeCollector::this_result& stepResult
        = result.get<typename StepVolumeCollector::result_type>();

    // Check the identified volumes
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (i > 0) {
        BOOST_TEST(stepResult.position[i].x() > 0.);
      }
      if (stepResult.position[i].x() >= 0.85 * units::_m
          && stepResult.position[i].x() < 0.95 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg2.name);
        BOOST_TEST(stepResult.volume[i]->material() == nullptr);
      } else if (stepResult.position[i].x() >= 1.05 * units::_m
                 && stepResult.position[i].x() < 1.15 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg1.name);
        BOOST_TEST(stepResult.volume[i]->material() != nullptr);
      } else if (stepResult.position[i].x() < 2. * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg.name);
        BOOST_TEST(stepResult.volume[i]->material() == nullptr);
      }
    }
  }

  BOOST_AUTO_TEST_CASE(BoxGeometryBuilderTest_confinedVolumes_edgecases)
  {
    // Production factory
    BoxGeometryBuilder bgb;

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg1;
    vCfg1.position = {1. * units::_m, 0., 0.};
    vCfg1.length   = {2. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg1.name     = "Test volume1";
    // Build and add 4 confined volumes
    // Volume that is missed and quite orthogonal to the starting position
    BoxGeometryBuilder::VolumeConfig cvCfg1;
    cvCfg1.position = {0.1 * units::_m, 0.4 * units::_m, 0.4 * units::_m};
    cvCfg1.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg1.name     = "Confined volume1";
    cvCfg1.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    // Volume that is missed but far away such that it may be hit
    BoxGeometryBuilder::VolumeConfig cvCfg2;
    cvCfg2.position = {1.9 * units::_m, -0.4 * units::_m, -0.4 * units::_m};
    cvCfg2.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg2.name     = "Confined volume2";
    // Volume that is hit but with identical boundary as its mother
    BoxGeometryBuilder::VolumeConfig cvCfg3;
    cvCfg3.position = {1.95 * units::_m, 0., 0.};
    cvCfg3.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg3.name     = "Confined volume3";
    // Volume to grind along the boundary
    BoxGeometryBuilder::VolumeConfig cvCfg4;
    cvCfg4.position = {1. * units::_m, 5. * units::_cm, 0.};
    cvCfg4.length   = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
    cvCfg4.name     = "Confined volume4";
    vCfg1.volumeCfg = {cvCfg1, cvCfg2, cvCfg3, cvCfg4};

    // Build a volume that confines another volume
    BoxGeometryBuilder::VolumeConfig vCfg2;
    vCfg2.position = {2.5 * units::_m, 0., 0.};
    vCfg2.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vCfg2.name     = "Test volume2";

    // Build detector
    BoxGeometryBuilder::Config config;
    config.position  = {1.5 * units::_m, 0., 0.};
    config.length    = {3. * units::_m, 1. * units::_m, 1. * units::_m};
    config.volumeCfg = {vCfg1, vCfg2};
    std::shared_ptr<TrackingGeometry> detector
        = bgb.buildTrackingGeometry(config);

    // Set propagator and navigator
    PropagatorOptions<ActionList<StepVolumeCollector>> propOpts;
    propOpts.maxStepSize = 10. * units::_mm;
    StraightLineStepper sls;
    Navigator           navi(detector);
    navi.resolvePassive   = true;
    navi.resolveMaterial  = true;
    navi.resolveSensitive = true;

    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        nullptr, startParams, startMom, 1.);

    // Launch and collect results
    const auto& result = prop.propagate(sbtp, propOpts);
    const StepVolumeCollector::this_result& stepResult
        = result.get<typename StepVolumeCollector::result_type>();

    // Check the identified volumes
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (i > 0) {
        BOOST_TEST(stepResult.position[i].x() > 0.);
      }
      if (stepResult.position[i].x() >= 0.95 * units::_m
          && stepResult.position[i].x() < 1.05 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg4.name);
      } else {
        if (stepResult.position[i].x() >= 1.9 * units::_m
            && stepResult.position[i].x() < 2. * units::_m) {
          BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg3.name);
        } else {
          if (stepResult.position[i].x() < 2. * units::_m) {
            BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg1.name);
          } else {
            if (stepResult.position[i].x() < 3. * units::_m)
              BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg2.name);
          }
        }
      }
    }
  }
}  // namespace Test
}  // namespace Acts
