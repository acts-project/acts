// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Extrapolator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "../Extrapolator/ExtrapolatorTestGeometry.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/Fitter/KalmanSequencer.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialCollector.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAbortConditions.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  using path_limit = detail::PathLimitReached;

  std::vector<std::unique_ptr<const Surface>> stepState;
  auto tGeometry = testGeometry<PlaneSurface>(stepState);

  using KalmanNavigator = Navigator<KalmanSequencer>;

  // get the navigator and provide the TrackingGeometry
  KalmanNavigator navigator(tGeometry);
  using BFieldType          = ConstantBField;
  using EigenStepperType    = EigenStepper<BFieldType>;
  using EigenPropagatorType = Propagator<EigenStepperType, KalmanNavigator>;

  const double        Bz = 2. * units::_T;
  BFieldType          bField(0, 0, Bz);
  EigenStepperType    estepper(bField);
  EigenPropagatorType epropagator(std::move(estepper), std::move(navigator));

  const int ntests    = 100;
  const int skip      = 0;
  bool      debugMode = false;

  // This test case checks Kalman fitter sequencing & reverse navigation
  BOOST_DATA_TEST_CASE(
      kalman_sequencer_test,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.25 * units::_GeV,
                                                        1.5 * units::_GeV)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
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

    if (index < skip) return;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    CurvilinearParameters start(nullptr, pos, mom, q);

    // Debug output
    using DebugOutput = detail::DebugOutputActor;
    using EndOfWorld  = detail::EndOfWorldReached;

    // Action list and abort list
    using ActionListType      = ActionList<DebugOutput>;
    using AbortConditionsType = AbortList<EndOfWorld>;

    PropagatorOptions<ActionListType, AbortConditionsType> options;
    options.maxStepSize = 25. * units::_cm;
    options.pathLimit   = 1500. * units::_mm;
    options.debug       = debugMode;

    const auto& status = epropagator.propagate(start, options);
    // get the backward output to the screen
    if (options.debug) {
      const auto& kalmanSequencerOutput
          = status.template get<DebugOutput::result_type>();
      std::cout << ">>> Output of the Kalman sequencer propagation "
                << std::endl;
      std::cout << kalmanSequencerOutput.debugString << std::endl;
    }

    BOOST_TEST(status.endParameters != nullptr);
    BOOST_TEST(status.endParameters->position().norm() < 0.25);
  }

}  // namespace Test
}  // namespace Acts
