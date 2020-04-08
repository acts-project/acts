// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <limits>
#include <random>
#include <vector>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Material/VolumeMaterialMapper.hpp"

using Acts::VectorHelpers::perp;

namespace Acts {

/// @brief create a small tracking geometry to map some dummy material on
std::shared_ptr<const TrackingGeometry> trackingGeometry() {
  using namespace Acts::UnitLiterals;

  BinUtility bu(5, 100, 1100, open, binR);
  bu += BinUtility(2, -M_PI, M_PI, closed, binPhi);
  bu += BinUtility(3, -1000, 1000, open, binZ);
  auto matProxy = std::make_shared<const ProtoVolumeMaterial>(bu);

  Logging::Level volumeLLevel = Logging::INFO;

  // tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      tvacConfig, getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));

  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config cvbConfig;
  cvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  cvbConfig.volumeName = "BeamPipe";
  cvbConfig.layerEnvelopeR = {100_mm, 1000_mm};
  cvbConfig.layerEnvelopeZ = 1000_mm;
  cvbConfig.buildToRadiusZero = false;
  cvbConfig.volumeMaterial = matProxy;
  cvbConfig.volumeSignature = 0;
  auto centralVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
      cvbConfig, getDefaultLogger("CentralVolumeBuilder", volumeLLevel));

  // create the bounds and the volume
  auto centralVolumeBounds =
      std::make_shared<const CylinderVolumeBounds>(100., 1000., 2000.);

  GeometryContext gCtx;
  auto centralVolume =
      centralVolumeBuilder->trackingVolume(gCtx, nullptr, centralVolumeBounds);

  return std::make_shared<const TrackingGeometry>(centralVolume);
}

std::shared_ptr<const TrackingGeometry> tGeometry = trackingGeometry();

namespace Test {

/// Test the filling and conversion
BOOST_AUTO_TEST_CASE(SurfaceMaterialMapper_tests) {
  /// We need a Navigator, Stepper to build a Propagator
  Navigator navigator(tGeometry);
  StraightLineStepper stepper;
  VolumeMaterialMapper::StraightLinePropagator propagator(std::move(stepper),
                                                          std::move(navigator));

  /// The config object
  Acts::VolumeMaterialMapper::Config vmmConfig;
  Acts::VolumeMaterialMapper vmMapper(vmmConfig, std::move(propagator));

  /// Create some contexts
  GeometryContext gCtx;
  MagneticFieldContext mfCtx;

  /// Now create the mapper state
  auto mState = vmMapper.createState(gCtx, mfCtx, *tGeometry);

  /// Test if this is not null
  BOOST_CHECK_EQUAL(mState.volumeMaterial.size(), 1u);
}

}  // namespace Test
}  // namespace Acts
