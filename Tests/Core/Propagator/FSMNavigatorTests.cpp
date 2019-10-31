// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Navigator Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <memory>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/FSMNavigator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

using VectorHelpers::perp;

namespace Test {
GeometryContext tgContext = GeometryContext();

/// This is a simple cache struct to mimic a Stepper
struct PseudoStepper {
  /// This is a simple cache struct to mimic the
  /// Stepper cache in the propagation
  struct State {
    /// Position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// Direction
    Vector3D dir = Vector3D(1., 0., 0.);

    /// Momentum
    double p;

    /// the navigation direction
    NavigationDirection navDir = forward;

    // accummulated path length cache
    double pathAccumulated = 0.;

    // adaptive sep size of the runge-kutta integration
    detail::ConstrainedStep stepSize =
        detail::ConstrainedStep(100 * units::_cm);
  };

  void step(State& sstate, double fraction) {
    // update the cache position
    double ssize = sstate.stepSize * fraction;
    std::cout << "PseudoStepper: Performing step with size: " << ssize
              << " along [" << sstate.dir.transpose() << "]: " << std::endl;
    std::cout << "               [" << sstate.pos.transpose();
    sstate.pos = sstate.pos + (sstate.stepSize * fraction) * sstate.dir;
    std::cout << "] -> [" << sstate.pos.transpose() << "]" << std::endl;
    // create navigation parameters
    return;
  }

  /// Global particle position accessor
  Vector3D position(const State& state) const { return state.pos; }

  /// Momentum direction accessor
  Vector3D direction(const State& state) const { return state.dir; }

  /// Momentum accessor
  Vector3D momentum(const State& state) const { return state.p * state.dir; }

  /// Charge access
  double charge(const State& /*state*/) const { return -1.; }

  /// Return a corrector
  static VoidIntersectionCorrector corrector(State& /*unused*/) {
    return VoidIntersectionCorrector();
  }

  bool surfaceReached(const State& state, const Surface* surface) const {
    return surface->isOnSurface(tgContext, position(state), direction(state),
                                true);
  }
};

/// This is a simple cache struct to mimic the
/// Propagator state
struct PropagatorState {
  /// emulate the options template
  struct Options {
    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;
  };

  /// Navigation cache: the start surface
  const Surface* startSurface = nullptr;

  /// Navigation cache: the current surface
  const Surface* currentSurface = nullptr;

  /// Navigation cache: the target surface
  const Surface* targetSurface = nullptr;
  bool targetReached = false;

  /// Give some options
  Options options;

  /// The Stepper state - internal statew of the Stepper
  PseudoStepper::State stepping;

  /// Navigation state - internal state of the Navigator
  FSMNavigator::state_type navigation;

  GeometryContext geoContext = GeometryContext();
};

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

BOOST_AUTO_TEST_CASE(Navigation) {
  FSMNavigator::Config cfg;
  cfg.trackingGeometry = tGeometry;
  cfg.resolveSensitive = true;
  cfg.resolveMaterial = true;
  cfg.resolvePassive = true;
  FSMNavigator navigator(cfg,
                         getDefaultLogger("FSMNavigator", Logging::VERBOSE));

  Vector3D position(0., 0., 0);
  Vector3D momentum(1., 1., 0);

  PseudoStepper stepper;
  PropagatorState state;

  auto step = [&](double fraction = 1.) {
    stepper.step(state.stepping, fraction);
    navigator.status(state, stepper);
    navigator.target(state, stepper);
  };

  state.stepping.pos = position;
  state.stepping.dir = momentum.normalized();

  navigator.status(state, stepper);

  // currentVolume has been set
  BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);

  BOOST_CHECK_EQUAL(state.navigation.currentVolume->volumeName(),
                    "BeamPipe::Barrel");

  // status is done, we should still be in initial state
  // propagator now calls target
  navigator.target(state, stepper);

  // do the first step
  step(0.5);
  step(1.0);
}

}  // namespace Test
}  // namespace Acts
