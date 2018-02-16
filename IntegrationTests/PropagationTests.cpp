// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE Propagator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank
#include <boost/test/data/test_case.hpp>

#include <cmath>

#include <boost/test/data/test_case.hpp>
#include "ACTS/ACTSVersion.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/Extrapolation/Wrapper.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/AtlasStepper.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/StrawSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "PropagationTestHelper.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using namespace propagation;

namespace IntegrationTest {

  typedef ConstantBField                BField_type;
  typedef EigenStepper<BField_type>     EigenStepper_type;
  typedef AtlasStepper<BField_type>     AtlasStepper_type;
  typedef Propagator<EigenStepper_type> EigenPropagator_type;
  typedef Propagator<AtlasStepper_type> AtlasPropagator_type;
  typedef RungeKuttaEngine<BField_type> PropagationEngine_type;

  typedef Wrapper<std::shared_ptr<PropagationEngine_type>>
      WrappedPropagator_type;

  // number of tests
  const int  ntests = 100000;
  const bool covtpr = true;

  // setup propagator with constant B-field
  const double         Bz = 2. * units::_T;
  BField_type          bField(0, 0, Bz);
  EigenStepper_type    estepper(bField);
  EigenPropagator_type epropagator(std::move(estepper));
  AtlasStepper_type    astepper(bField);
  AtlasPropagator_type apropagator(std::move(astepper));
  auto                 bFieldPtr = std::make_shared<const BField_type>(bField);
  auto                 wConfig   = PropagationEngine_type::Config(bFieldPtr);
  auto                 wegine    = std::make_shared<PropagationEngine_type>(
      wConfig,
      Acts::getDefaultLogger("RungeKuttaEngine", Acts::Logging::INFO));
  WrappedPropagator_type wpropagator(wegine);

  // The constant field test
  /// test forward propagation in constant magnetic field
  BOOST_DATA_TEST_CASE(
      constant_bfieldorward_propagation_,
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
    double dcharge = -1 + 2 * charge;
    // constant field propagation atlas stepper
    auto aposition = constant_field_propagation(
        apropagator, pT, phi, theta, dcharge, index, Bz);
    // constant field propagation eigen stepper
    auto eposition = constant_field_propagation(
        epropagator, pT, phi, theta, dcharge, index, Bz);
    // constant field runge kutta engine - not yet at same accuracy
    auto wposition = constant_field_propagation(
        wpropagator, pT, phi, theta, dcharge, index, Bz, 10. * units::_um);
    // check consistency
    BOOST_CHECK(eposition.isApprox(aposition));
    BOOST_CHECK(eposition.isApprox(wposition, 1e-3));
  }

  // The actual test - needs to be included to avoid
  // template inside template definition through boost
  // #include "PropagationTestBase.hpp"

}  // namespace Test

}  // namespace Acts
