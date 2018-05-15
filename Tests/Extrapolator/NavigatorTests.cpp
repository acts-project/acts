// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE MaterialCollection Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include <memory>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolator/Navigator.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// This is a simple cache struct to mimic the
  /// Propagator cache
  struct PropagatorCache
  {
    /// emulate the options template
    struct Options
    {
      /// Debug output
      /// the string where debug messages are stored (optionally)
      bool        debug       = false;
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
    bool           targetReached = false;

    /// Give some options
    Options options;
  };

  /// This is a simple cache struct to mimic the
  /// Stepper cache in the propagation
  struct StepperCache
  {

    /// Access method to satisify TrackingVolume interface
    const Vector3D&
    position() const
    {
      return pos;
    }

    /// Access method to satisify TrackingVolume interface
    const Vector3D&
    momentum() const
    {
      return dir;
    }

    /// Position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// and mumentum
    Vector3D dir = Vector3D(1., 0., 0.);

    /// the navigation direction
    NavigationDirection navDir = forward;

    // accummulated path length cache
    double accumulatedPath = 0.;

    // adaptive sep size of the runge-kutta integration
    cstep stepSize = 100 * units::_cm;
  };

  template <typename stepper_state_t>
  NavigationParameters
  step(stepper_state_t& stepState)
  {
    // update the cache position
    stepState.pos = stepState.pos + stepState.stepSize * stepState.dir;
    // create navigation parameters
    return NavigationParameters(stepState.pos, stepState.dir);
  }

  // the surface cache & the creation of the geometry
  std::vector<std::unique_ptr<const Surface>> surfaceCache;
  auto tGeometry = testGeometry<PlaneSurface>(surfaceCache);

  // the debug boolean
  bool debug = true;

  BOOST_AUTO_TEST_CASE(Navigator_methods)
  {

    // create a navigator
    Navigator navigator;
    navigator.trackingGeometry = tGeometry;
    navigator.debug            = true;
    navigator.collectSensitive = true;
    navigator.collectMaterial  = true;
    navigator.collectPassive   = false;

    // position and direction vector
    Vector3D position(0., 0., 0);
    Vector3D momentum(1., 1., 0);

    // Navigation parameter
    NavigationParameters navPar(position, momentum);

    // the result and the cache
    Navigator::result_type result;

    // the propagator cache
    PropagatorCache propState;
    propState.options.debug = debug;

    // the stepper cache
    StepperCache stepState;
    stepState.pos = position;
    stepState.dir = momentum.unit();

    // (1) Initialization navigation from start point
    // - this will call resolveLayers() as well
    // - it will not return to the stepper
    BOOST_TEST(navigator.initialize(navPar, propState, stepState, result));
    // check that the currentVolume is set
    BOOST_TEST((result.currentVolume != nullptr));
    // check that the currentVolume is the startVolume
    BOOST_TEST((result.currentVolume == result.startVolume));
    // check that the currentSurface is reset to:
    BOOST_TEST((propState.currentSurface == nullptr));
    // one layer has been found
    BOOST_TEST((result.navLayers.size() == 1));
    // the iterator should point to it
    BOOST_TEST((result.navLayerIter == result.navLayers.begin()));
    // cahce the beam pipe radius
    double beamPipeRadius = result.navLayerIter->intersection.position.perp();
    // step size has been updated
    BOOST_TEST((stepState.stepSize == beamPipeRadius));
    if (debug) {
      std::cout << "<<< Test 1a >>> initialize at (0,0,0) " << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar              = step(stepState);
    Vector3D onBeamPipe = navPar.position();

    // (2) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        !(navigator.handleSurfaces(navPar, propState, stepState, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handleLayers(navPar, propState, stepState, result)));
    BOOST_TEST((propState.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(
        navigator.handleBoundaries(navPar, propState, stepState, result));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1b >>> handleLayers, and set step to boundary  "
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar              = step(stepState);
    Vector3D onBoundary = navPar.position();

    // (3) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        !(navigator.handleSurfaces(navPar, propState, stepState, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handleLayers(navPar, propState, stepState, result)));
    BOOST_TEST((propState.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(
        navigator.handleBoundaries(navPar, propState, stepState, result));

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1c >>> advance to boundary, initialize layers."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar                 = step(stepState);
    Vector3D on1stApproach = navPar.position();

    // (4) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        !(navigator.handleSurfaces(navPar, propState, stepState, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, propState, stepState, result)));

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1d >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar            = step(stepState);
    Vector3D on1stSf1 = navPar.position();

    // (5) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1e >>> advance to surface, update." << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (6) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1f >>> advance to next surface, update."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (7) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1g >>> advance to next surface, update."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (8) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1h >>> advance to next surface, update."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (9) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        !(navigator.handleSurfaces(navPar, propState, stepState, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1i >>> advance to last surface, switch layer."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (10) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        !(navigator.handleSurfaces(navPar, propState, stepState, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, propState, stepState, result)));

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1j >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (11) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1k >>> advance to surface, update." << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // positive return: do the step
    navPar = step(stepState);

    // (12) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // handleSurfaces should return false
    BOOST_TEST(
        (navigator.handleSurfaces(navPar, propState, stepState, result)));
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1l >>> advance to next surface, update."
                << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // Initialize navigation from beam pipe
    result = Navigator::result_type();
    // recreate cache
    stepState     = StepperCache();
    stepState.pos = onBeamPipe;
    stepState.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(stepState.pos, momentum);
    BOOST_TEST(navPar.position() == onBeamPipe);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // step size has been updated
    if (debug) {
      std::cout << "<<< Test 2 >>> initialize at BeamPipe " << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // Initialize navigation from boundary surface
    result = Navigator::result_type();
    // recreate cache
    stepState     = StepperCache();
    stepState.pos = onBoundary;
    stepState.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(stepState.pos, momentum);
    BOOST_TEST(navPar.position() == onBoundary);
    BOOST_TEST((navigator.initialize(navPar, propState, stepState, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 3 >>> initialize from Boundary " << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // Initialize navigation from approach surface of first layer
    result = Navigator::result_type();
    // recreate cache
    stepState     = StepperCache();
    stepState.pos = on1stApproach;
    stepState.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(stepState.pos, momentum);
    BOOST_TEST(navPar.position() == on1stApproach);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 4 >>> initialize from Approach " << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }

    // Initialize navigation from 1st surface of
    result = Navigator::result_type();
    // recreate cache
    stepState     = StepperCache();
    stepState.pos = on1stSf1;
    stepState.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(stepState.pos, momentum);
    BOOST_TEST(navPar.position() == on1stSf1);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, propState, stepState, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 5 >>> initialize from 1st Surface " << std::endl;
      std::cout << propState.options.debugString << std::endl;
      propState.options.debugString = "";
    }
  }
}

}  // end of namespace Acts
