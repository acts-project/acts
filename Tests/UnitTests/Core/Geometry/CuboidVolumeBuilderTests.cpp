// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <vector>

#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

struct StepVolumeCollector {
  ///
  /// @brief Data container for result analysis
  ///
  struct this_result {
    // Position of the propagator after each step
    std::vector<Vector3D> position;
    // Volume of the propagator after each step
    std::vector<TrackingVolume const*> volume;
  };

  using result_type = this_result;

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    result.position.push_back(stepper.position(state.stepping));
    result.volume.push_back(state.navigation.currentVolume);
  }
};

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
    MaterialProperties matProp(makeBeryllium(), 0.5_mm);
    cfg.surMat = std::make_shared<HomogeneousSurfaceMaterial>(matProp);

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
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());

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
    MaterialProperties matProp(makeBeryllium(), 0.5_mm);
    cfg.surMat = std::make_shared<HomogeneousSurfaceMaterial>(matProp);

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

/*
BOOST_AUTO_TEST_CASE(CuboidVolumeBuilderTest_confinedVolumes) {
  // Production factory
  CuboidVolumeBuilder cvb;

  // Create a test context
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();

  // Build a volume that confines another volume
  CuboidVolumeBuilder::VolumeConfig vCfg;
  vCfg.position = {1. * units::_m, 0., 0.};
  vCfg.length = {2. * units::_m, 1. * units::_m, 1. * units::_m};
  vCfg.name = "Test volume";
  // Build and add 2 confined volumes
  CuboidVolumeBuilder::VolumeConfig cvCfg1;
  cvCfg1.position = {1.1 * units::_m, 0., 0.};
  cvCfg1.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg1.name = "Confined volume1";
  cvCfg1.volumeMaterial =
      std::shared_ptr<const IVolumeMaterial>(new HomogeneousVolumeMaterial(
          Material(352.8, 407., 9.012, 4., 1.848e-3)));
  CuboidVolumeBuilder::VolumeConfig cvCfg2;
  cvCfg2.position = {0.9 * units::_m, 0., 0.};
  cvCfg2.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg2.name = "Confined volume2";
  vCfg.volumeCfg = {cvCfg1, cvCfg2};

  // Build detector
  CuboidVolumeBuilder::Config config;
  config.position = {1. * units::_m, 0., 0.};
  config.length = {2. * units::_m, 1. * units::_m, 1. * units::_m};
  config.volumeCfg = {vCfg};

  cvb.setConfig(config);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  // Test that the right volume is selected
  BOOST_TEST(detector->lowestTrackingVolume(tgContext, {1. * units::_m, 0., 0.})
                 ->volumeName() == vCfg.name);
  BOOST_TEST(
      detector->lowestTrackingVolume(tgContext, {1.1 * units::_m, 0., 0.})
          ->volumeName() == cvCfg1.name);
  BOOST_TEST(
      detector->lowestTrackingVolume(tgContext, {0.9 * units::_m, 0., 0.})
          ->volumeName() == cvCfg2.name);

  // Set propagator and navigator
  PropagatorOptions<ActionList<StepVolumeCollector>> propOpts(tgContext,
                                                              mfContext);
  propOpts.maxStepSize = 10. * units::_mm;
  StraightLineStepper sls;
  Navigator navi(detector);
  navi.resolvePassive = true;
  navi.resolveMaterial = true;
  navi.resolveSensitive = true;

  Propagator<StraightLineStepper, Navigator> prop(sls, navi);

  // Set initial parameters for the particle track
  Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
      std::nullopt, startParams, startMom, 1., 0.);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepVolumeCollector::this_result& stepResult =
      result.get<typename StepVolumeCollector::result_type>();

  // Check the identified volumes
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    if (i > 0) {
      BOOST_TEST(stepResult.position[i].x() > 0.);
    }
    if (stepResult.position[i].x() >= 0.85 * units::_m &&
        stepResult.position[i].x() < 0.95 * units::_m) {
      BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg2.name);
      BOOST_TEST(stepResult.volume[i]->volumeMaterial() == nullptr);
    } else {
      if (stepResult.position[i].x() >= 1.05 * units::_m &&
          stepResult.position[i].x() < 1.15 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg1.name);
        BOOST_TEST(stepResult.volume[i]->volumeMaterial() != nullptr);
      } else {
        if (stepResult.position[i].x() < 2. * units::_m) {
          BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg.name);
          BOOST_TEST(stepResult.volume[i]->volumeMaterial() == nullptr);
        }
      }
    }
  }
}
*/
/** As discussed with the author, disabled for the moment
BOOST_AUTO_TEST_CASE(CuboidVolumeBuilderTest_confinedVolumes_edgecases) {
  // Production factory
  CuboidVolumeBuilder cvb;

  // Create a test context
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();

  // Build a volume that confines another volume
  CuboidVolumeBuilder::VolumeConfig vCfg1;
  vCfg1.position = {1. * units::_m, 0., 0.};
  vCfg1.length = {2. * units::_m, 1. * units::_m, 1. * units::_m};
  vCfg1.name = "Test volume1";
  // Build and add 4 confined volumes
  // Volume that is missed and quite orthogonal to the starting position
  CuboidVolumeBuilder::VolumeConfig cvCfg1;
  cvCfg1.position = {0.1 * units::_m, 0.4 * units::_m, 0.4 * units::_m};
  cvCfg1.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg1.name = "Confined volume1";
  cvCfg1.volumeMaterial =
      std::shared_ptr<const IVolumeMaterial>(new HomogeneousVolumeMaterial(
          Material(352.8, 407., 9.012, 4., 1.848e-3)));
  // Volume that is missed but far away such that it may be hit
  CuboidVolumeBuilder::VolumeConfig cvCfg2;
  cvCfg2.position = {1.9 * units::_m, -0.4 * units::_m, -0.4 * units::_m};
  cvCfg2.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg2.name = "Confined volume2";
  // Volume that is hit but with identical boundary as its mother
  // TODO: Moved slightly inside the volume since otherwise the navigation
  // breaks due to overlapping boundary surfaces. The corresponding test below
  // is changed accordingly.
  CuboidVolumeBuilder::VolumeConfig cvCfg3;
  cvCfg3.position = {1.9 * units::_m, 0., 0.};
  cvCfg3.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg3.name = "Confined volume3";
  // Volume to grind along the boundary
  CuboidVolumeBuilder::VolumeConfig cvCfg4;
  cvCfg4.position = {1. * units::_m, 5. * units::_cm, 0.};
  cvCfg4.length = {10. * units::_cm, 10. * units::_cm, 10. * units::_cm};
  cvCfg4.name = "Confined volume4";
  vCfg1.volumeCfg = {cvCfg1, cvCfg2, cvCfg3, cvCfg4};

  // Build a volume that confines another volume
  CuboidVolumeBuilder::VolumeConfig vCfg2;
  vCfg2.position = {2.5 * units::_m, 0., 0.};
  vCfg2.length = {1. * units::_m, 1. * units::_m, 1. * units::_m};
  vCfg2.name = "Test volume2";

  // Build detector
  CuboidVolumeBuilder::Config config;
  config.position = {1.5 * units::_m, 0., 0.};
  config.length = {3. * units::_m, 1. * units::_m, 1. * units::_m};
  config.volumeCfg = {vCfg1, vCfg2};

  cvb.setConfig(config);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  using DebugOutput = Acts::DebugOutputActor;

  // Set propagator and navigator
  PropagatorOptions<ActionList<StepVolumeCollector,DebugOutput>>
propOpts(tgContext, mfContext);

  propOpts.debug = true;
  propOpts.maxStepSize = 10. * units::_mm;
  StraightLineStepper sls;

  Navigator navi(detector);
  navi.resolvePassive = true;
  navi.resolveMaterial = true;
  navi.resolveSensitive = true;

  Propagator<StraightLineStepper, Navigator> prop(sls, navi);

  // Set initial parameters for the particle track
  Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
      std::nullopt, startParams, startMom, 1., 0.);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepVolumeCollector::this_result& stepResult =
      result.get<typename StepVolumeCollector::result_type>();

  // Screen output
  if (propOpts.debug) {
    const auto debugString =
        result.template get<DebugOutput::result_type>().debugString;
    std::cout << debugString << std::endl;
  }

  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    // Check the movement in the right direction
    if (i > 0) {
      BOOST_TEST(stepResult.position[i].x() > 0.);
    }
    // Check the identified volumes
    if (stepResult.position[i].x() >= 0.95 * units::_m &&
        stepResult.position[i].x() < 1.05 * units::_m) {
      BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg4.name);
    } else {
      if (stepResult.position[i].x() >= 1.85 * units::_m &&
          stepResult.position[i].x() < 1.95 * units::_m) {
        BOOST_TEST(stepResult.volume[i]->volumeName() == cvCfg3.name);
      } else {
        if (stepResult.position[i].x() < 2. * units::_m) {
          BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg1.name);
        } else {
          if (stepResult.position[i].x() < 3. * units::_m) {
            BOOST_TEST(stepResult.volume[i]->volumeName() == vCfg2.name);
          }
        }
      }
    }
  }
}
*/

}  // namespace Test
}  // namespace Acts
