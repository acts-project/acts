// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Propagator Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <cmath>
#include <iostream>

#include "Acts/ActsVersion.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "PropagationTestHelper.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace IntegrationTest {

  using BFieldType       = ConstantBField;
  using EigenStepperType = EigenStepper<BFieldType>;
  using DenseStepperType
      = EigenStepper<BFieldType,
                     VoidIntersectionCorrector,
                     StepperExtensionList<DenseEnvironmentExtension>>;
  using AtlasStepperType    = AtlasStepper<BFieldType>;
  using EigenPropagatorType = Propagator<EigenStepperType>;
  using DensePropagatorType = Propagator<DenseStepperType, Navigator>;
  using AtlasPropagatorType = Propagator<AtlasStepperType>;

  // number of tests
  const int  ntests = 100;
  const int  skip   = 0;
  const bool covtpr = true;
  const bool debug  = false;

  // setup propagator with constant B-field
  const double        Bz = 2. * units::_T;
  BFieldType          bField(0, 0, Bz);
  EigenStepperType    estepper(bField);
  DenseStepperType    dstepper(bField);
  EigenPropagatorType epropagator(std::move(estepper));
  AtlasStepperType    astepper(bField);
  AtlasPropagatorType apropagator(std::move(astepper));

  DensePropagatorType
  setupDensePropagator()
  {
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {1.5 * units::_m, 0., 0.};
    vConf.length   = {3. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf.material = std::make_shared<const Material>(
        Material(352.8, 407., 9.012, 4., 1.848e-3));
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {1.5 * units::_m, 0., 0.};
    conf.length   = {3. * units::_m, 1. * units::_m, 1. * units::_m};
    CuboidVolumeBuilder             cvb(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> detector = tgb.trackingGeometry();
    Navigator                               navi(detector);
    return DensePropagatorType(dstepper, std::move(navi));
  }

  // The constant field test
  /// test forward propagation in constant magnetic field
  BOOST_DATA_TEST_CASE(
      constant_bfieldforward_propagation_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    if (index < skip) {
      return;
    }

    double dcharge = -1 + 2 * charge;
    // constant field propagation atlas stepper
    auto aposition = constant_field_propagation(
        apropagator, pT, phi, theta, dcharge, index, Bz);
    // constant field propagation eigen stepper
    auto eposition = constant_field_propagation(
        epropagator, pT, phi, theta, dcharge, index, Bz);
    // check consistency
    CHECK_CLOSE_REL(eposition, aposition, 1e-6);
  }

// The actual test - needs to be included to avoid
// template inside template definition through boost
#include "PropagationTestBase.hpp"

}  // namespace Test

}  // namespace Acts
