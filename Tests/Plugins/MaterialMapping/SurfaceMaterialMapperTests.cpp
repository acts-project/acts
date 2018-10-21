// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE SurfaceMaterialMapper Tests
#include <boost/test/included/unit_test.hpp>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterialProxy.hpp"
#include "Acts/Plugins/MaterialMapping/SurfaceMaterialMapper.hpp"
#include "Acts/Tools/CylinderVolumeBuilder.hpp"
#include "Acts/Tools/CylinderVolumeHelper.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Tools/PassiveLayerBuilder.hpp"
#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

namespace Acts {

/// @brief create a small tracking geometry to map some dummy material on
std::shared_ptr<const TrackingGeometry>
trackingGeometry()
{

  BinUtility zbinned(8, -40, 40, open, binZ);
  auto       matProxy = std::make_shared<const SurfaceMaterialProxy>(zbinned);

  Logging::Level surfaceLLevel = Logging::INFO;
  Logging::Level layerLLevel   = Logging::INFO;
  Logging::Level volumeLLevel  = Logging::INFO;

  // configure surface array creator
  auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
      getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator            = std::make_shared<const LayerCreator>(
      lcConfig, getDefaultLogger("LayerCreator", layerLLevel));
  // configure the layer array creator
  auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
      getDefaultLogger("LayerArrayCreator", layerLLevel));

  // tracking volume array creator
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator          = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

  PassiveLayerBuilder::Config layerBuilderConfig;
  layerBuilderConfig.layerIdentification     = "CentralBarrel";
  layerBuilderConfig.centralLayerRadii       = {10., 20., 30.};
  layerBuilderConfig.centralLayerHalflengthZ = {40., 40., 40.};
  layerBuilderConfig.centralLayerThickness   = {1., 1., 1.};
  layerBuilderConfig.centralLayerMaterial    = {matProxy, matProxy, matProxy};
  auto layerBuilder = std::make_shared<const PassiveLayerBuilder>(
      layerBuilderConfig,
      getDefaultLogger("CentralBarrelBuilder", layerLLevel));
  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config cvbConfig;
  cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  cvbConfig.volumeName           = "BeamPipe";
  cvbConfig.layerBuilder         = layerBuilder;
  cvbConfig.layerEnvelopeR       = {1. * units::_mm, 1. * units::_mm};
  cvbConfig.buildToRadiusZero    = true;
  cvbConfig.volumeSignature      = 0;
  auto centralVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
      cvbConfig, getDefaultLogger("CentralVolumeBuilder", volumeLLevel));

  // create the bounds and the volume
  auto centralVolumeBounds
      = std::make_shared<const CylinderVolumeBounds>(0., 40., 110.);

  auto centralVolume
      = centralVolumeBuilder->trackingVolume(nullptr, centralVolumeBounds);

  return std::make_shared<const TrackingGeometry>(centralVolume);
}

std::shared_ptr<const TrackingGeometry> tGeometry = trackingGeometry();

namespace Test {

  /// Test the filling and conversion
  BOOST_AUTO_TEST_CASE(SurfaceMaterialMapper_tests)
  {

    /// We need a Navigator, Stepper to build a Propagator
    Navigator                                     navigator(tGeometry);
    StraightLineStepper                           stepper;
    SurfaceMaterialMapper::StraightLinePropagator propagator(
        std::move(stepper), std::move(navigator));

    /// The config object
    SurfaceMaterialMapper::Config smmConfig;
    SurfaceMaterialMapper         smMapper(smmConfig, std::move(propagator));

    /// Now create the mapper state
    auto mState = smMapper.createState(*tGeometry);

    /// Test if this is not null
    BOOST_CHECK_EQUAL(mState.accumulatedMaterial.size(), 3);

    // material properties
    MaterialProperties a(1., 1., 1., 1., 1., 1.);
    // and vacuum
    MaterialProperties v(1.);

    // we shoot under an angle of
    double cotan_theta_03_13_24 = 1.25 / 3.;
    // path scaled material
    MaterialProperties a_theta_03_13_24(a);
    a *= 1. / sin(atan2(3, 1.25));
    RecordedMaterialProperties m03{a, Vector3D(1., 0., cotan_theta_03_13_24)};
    RecordedMaterialProperties m13{a,
                                   Vector3D(2., 0., 2 * cotan_theta_03_13_24)};
    RecordedMaterialProperties m24{a,
                                   Vector3D(2., 0., 2 * cotan_theta_03_13_24)};
    std::vector<RecordedMaterialProperties> rmps = {m03, m13, m24};

    RecordedMaterialTrack rmt031324(
        Vector3D(0., 0., 0.), m03.second.normalized(), rmps);

    smMapper.mapMaterialTrack(mState, rmt031324);

    RecordedMaterialProperties m02{a, Vector3D(1., 0., -cotan_theta_03_13_24)};
    RecordedMaterialProperties m21{a,
                                   Vector3D(2., 0., -2 * cotan_theta_03_13_24)};

    RecordedMaterialTrack rmt02xx21(
        Vector3D(0., 0., 0.), m02.second.normalized(), rmps);

    smMapper.mapMaterialTrack(mState, rmt02xx21);
  }

}  // namespace Test

}  // namespace Acts
