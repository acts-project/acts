// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE SimpleGeometryTest

#include <boost/test/included/unit_test.hpp>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Tools/CylinderVolumeBuilder.hpp"
#include "Acts/Tools/CylinderVolumeHelper.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Tools/PassiveLayerBuilder.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  /// @brief Unit test for a three layer detector parameters
  /// Testing the Tool chain in the geometry building process
  ///
  BOOST_AUTO_TEST_CASE(SimpleGeometryTest)
  {

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
    auto tVolumeArrayCreator
        = std::make_shared<const TrackingVolumeArrayCreator>(
            getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
    // configure the cylinder volume helper
    CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator          = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
        cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

    // ----------------- build a beam pipe -----------------------------------
    PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii   = std::vector<double>(1, 3. * units::_mm);
    bplConfig.centralLayerHalflengthZ
        = std::vector<double>(1, 40. * units::_mm);
    bplConfig.centralLayerThickness = std::vector<double>(1, 0.8 * units::_mm);
    auto beamPipeBuilder = std::make_shared<const PassiveLayerBuilder>(
        bplConfig, getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
    // create the volume for the beam pipe
    CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName           = "BeamPipe";
    bpvConfig.layerBuilder         = beamPipeBuilder;
    bpvConfig.layerEnvelopeR       = {1. * units::_mm, 1. * units::_mm};
    bpvConfig.buildToRadiusZero    = true;
    bpvConfig.volumeSignature      = 0;
    auto beamPipeVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
        bpvConfig, getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));

    PassiveLayerBuilder::Config layerBuilderConfig;
    layerBuilderConfig.layerIdentification = "CentralBarrel";
    layerBuilderConfig.centralLayerRadii
        = {10. * units::_mm, 20. * units::_mm, 30. * units::_mm};
    layerBuilderConfig.centralLayerHalflengthZ
        = {40. * units::_mm, 40. * units::_mm, 40. * units::_mm};
    layerBuilderConfig.centralLayerThickness
        = {1. * units::_mm, 1. * units::_mm, 1. * units::_mm};
    auto layerBuilder = std::make_shared<const PassiveLayerBuilder>(
        layerBuilderConfig,
        getDefaultLogger("CentralBarrelBuilder", layerLLevel));
    // create the volume for the central barrel
    CylinderVolumeBuilder::Config cvbConfig;
    cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    cvbConfig.volumeName           = "CentralBarrel";
    cvbConfig.layerBuilder         = layerBuilder;
    cvbConfig.layerEnvelopeR       = {1. * units::_mm, 1. * units::_mm};
    cvbConfig.buildToRadiusZero    = false;
    cvbConfig.volumeSignature      = 0;
    auto centralVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
        cvbConfig, getDefaultLogger("CentralVolumeBuilder", volumeLLevel));

    // Make the TrackingGeometry Builder
    TrackingGeometryBuilder::Config tgbConfig;
    tgbConfig.trackingVolumeBuilders
        = {beamPipeVolumeBuilder, centralVolumeBuilder};
    tgbConfig.trackingVolumeHelper = cylinderVolumeHelper;

    TrackingGeometryBuilder tgBuilder(tgbConfig);
    auto                    tGeometry = tgBuilder.trackingGeometry();

    BOOST_CHECK(tGeometry != nullptr);
  }
}  // namespace Test
}  // namespace Acts
