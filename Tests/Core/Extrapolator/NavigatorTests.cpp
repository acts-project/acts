// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Navigator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include <memory>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  /// This is a simple cache struct to mimic the
  /// Propagator cache
  struct PropagatorState
  {

    /// This is a simple cache struct to mimic the
    /// Stepper cache in the propagation
    struct StepperState
    {

      /// Access method to position
      const Vector3D&
      position() const
      {
        return pos;
      }

      /// Access method to momentum
      const Vector3D&
      momentum() const
      {
        return dir;
      }

      /// Access method to direction
      const Vector3D&
      direction() const
      {
        return dir;
      }

      /// Return a corrector
      VoidCorrector
      corrector() const
      {
        return VoidCorrector();
      }

      /// Position
      Vector3D pos = Vector3D(0., 0., 0.);

      /// and mumentum
      Vector3D dir = Vector3D(1., 0., 0.);

      /// the navigation direction
      NavigationDirection navDir = forward;

      // accummulated path length cache
      double pathAccumulated = 0.;

      // adaptive sep size of the runge-kutta integration
      cstep stepSize = cstep(100 * units::_cm);
    };

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

    /// The Stepper state - internal statew of the Stepper
    StepperState stepping;

    /// Navigation state - internal state of the Navigator
    Navigator::state_type navigation;
  };

  template <typename stepper_state_t>
  void
  step(stepper_state_t& sstate)
  {
    // update the cache position
    sstate.pos = sstate.pos + sstate.stepSize * sstate.dir;
    // create navigation parameters
    return;
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
    navigator.initialStepFactor = 1.;
    navigator.trackingGeometry  = tGeometry;
    navigator.resolveSensitive  = true;
    navigator.resolveMaterial   = true;
    navigator.resolvePassive    = false;

    // position and direction vector
    Vector3D position(0., 0., 0);
    Vector3D momentum(1., 1., 0);

    // the propagator cache
    PropagatorState state;
    state.options.debug = debug;

    // the stepper cache
    state.stepping.pos = position;
    state.stepping.dir = momentum.unit();

    auto navCorr = state.stepping.corrector();

    // foward navigation ----------------------------------------------
    if (debug) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<< FORWARD NAVIGATION >>>>>>>>>>>>>>>>>>"
                << std::endl;
    }

    // (1) Initialization navigation from start point
    // - this will call resolveLayers() as well
    // - and thus should call a return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == true);
    // check that the currentVolume is set
    BOOST_TEST((state.navigation.currentVolume != nullptr));
    // check that the currentVolume is the startVolume
    BOOST_TEST(
        (state.navigation.currentVolume == state.navigation.startVolume));
    // check that the currentSurface is reset to:
    BOOST_TEST((state.navigation.currentSurface == nullptr));
    // one layer has been found
    BOOST_TEST((state.navigation.navLayers.size() == 1));
    // the iterator should point to it
    BOOST_TEST(
        (state.navigation.navLayerIter == state.navigation.navLayers.begin()));
    // cache the beam pipe radius
    double beamPipeRadius
        = state.navigation.navLayerIter->intersection.position.perp();
    // step size has been updated
    BOOST_TEST((state.stepping.stepSize == beamPipeRadius));
    if (debug) {
      std::cout << "<<< Test 1a >>> initialize at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);
    Vector3D onBeamPipe = state.stepping.position();

    // (2) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : do not return to the stepper, no surfaces set
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // handle_layer should store the layer surface, but return false
    BOOST_TEST(navigator.handleLayers(state, navCorr) == false);
    BOOST_TEST((state.navigation.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(navigator.handleBoundaries(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1b >>> handleLayers, and set step to boundary at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);
    Vector3D onBoundary = state.stepping.position();

    // (3) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : do not return to the stepper
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(navigator.handleLayers(state, navCorr) == false);
    BOOST_TEST((state.navigation.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(navigator.handleBoundaries(state, navCorr) == true);

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1c >>> advance to boundary, initialize layers at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);
    Vector3D on1stApproach = state.stepping.position();

    // (4) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : do not return to the stepper
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // handleLayers should sore the layer surface, return to stepper
    BOOST_TEST(navigator.handleLayers(state, navCorr) == true);

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1d >>> advance to layer, initialize surfaces at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);
    Vector3D on1stSf1 = state.stepping.position();

    // (5) re-entering navigator:
    // initialize: do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1e >>> advance to surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (6) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) = true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1f >>> advance to next surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (7) re-entering navigator:
    // initialize: do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1g >>> advance to next surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (8) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1h >>> advance to next surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (9) re-entering navigator:
    // initialize: do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : no more surfaces, do not return to the stepper
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // handleLayers : set step size towards it and return
    BOOST_TEST(navigator.handleLayers(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1i >>> advance to last surface, switch layer at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (10) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : no surfaces, do not return to the stepper
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(navigator.handleLayers(state, navCorr) == true);

    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1j >>> advance to layer, initialize surfaces at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (11) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == true);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1k >>> advance to surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // (12) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // handleSurfaces : set step size and return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == true);
    // check that we have the current surface
    BOOST_TEST(state.navigation.currentSurface);
    // the iterator should point to it
    if (debug) {
      std::cout << "<<< Test 1l >>> advance to next surface, update at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // remember the end parameters for the backward navigation
    Vector3D eposition = state.stepping.position();
    // remember the end surface for the backward navidation
    const Surface* esurface = state.navigation.currentSurface;

    // backward navigation ----------------------------------------------
    if (debug) {
      std::cout
          << "<<<<<<<<<<<<<<<<<<<<< BACKWARD NAVIGATION >>>>>>>>>>>>>>>>>>"
          << std::endl;
    }

    // re-initialize the propagator state
    state                         = PropagatorState();
    state.stepping.navDir         = backward;
    state.stepping.stepSize       = cstep(-100 * units::_cm);
    state.stepping.dir            = momentum.unit();
    state.stepping.pos            = eposition;
    state.options.debug           = debug;
    state.navigation.startSurface = esurface;

    // initialize the navigator
    // this one avoids the stepping towards layer
    // (1) Initialization navigation from start point
    // - this will call resolveLayers() as well, which for a
    //   start layer will call 'resolveSurfaces'
    // hence return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == true);

    // check that the currentVolume is set
    BOOST_TEST(state.navigation.currentVolume);
    // check that the currentVolume is the startVolume
    BOOST_TEST(
        (state.navigation.currentVolume == state.navigation.startVolume));
    // one layer has been found
    BOOST_TEST((state.navigation.navLayers.size()));
    if (debug) {
      std::cout << "<<< Test -1a >>> initialize at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
    // positive return: do the step
    step(state.stepping);

    // (3) re-entering navigator:
    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // there are no more surfaces after this one, do not return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    BOOST_TEST(navigator.handleLayers(state, navCorr) == true);
    if (debug) {
      std::cout << "<<< Test -1b >>> handle layer at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    // initialize : do not return to the stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // there are no more surfaces after this one, do not return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    BOOST_TEST(navigator.handleLayers(state, navCorr) == true);
    if (debug) {
      std::cout << "<<< Test -1e >>> handle layer at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step
    step(state.stepping);

    std::vector<std::string> ssteps = {"f", "g", "h", "i"};

    // go backwards throught the surfaces
    for (auto& s : ssteps) {
      BOOST_TEST(navigator.initialize(state, navCorr) == false);
      // handleSurfaces : set step size and return
      bool returnToStepper = navigator.handleSurfaces(state, navCorr);
      BOOST_TEST(returnToStepper == (true && s != "i"));
      // the iterator should point to it
      if (debug) {
        std::cout << "<<< Test -1";
        std::cout << s;
        std::cout << " >>> advance to next surface, update at "
                  << toString(state.stepping.position()) << std::endl;
        std::cout << state.options.debugString << std::endl;
        state.options.debugString = "";
      }

      // positive return: do the step
      if (returnToStepper) {
        step(state.stepping);
      }
    }

    // handle layers should give flag to move on to boundaries
    BOOST_TEST(navigator.handleLayers(state, navCorr) == false);
    BOOST_TEST(navigator.handleBoundaries(state, navCorr) == true);
    if (debug) {
      std::cout << "<<< Test -1j >>> advance to boundary surface from "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
    // positive return: do the step
    step(state.stepping);

    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // there are no more surfaces after this one, do not return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    BOOST_TEST(navigator.handleLayers(state, navCorr) == false);
    BOOST_TEST(navigator.handleBoundaries(state, navCorr) == true);
    if (debug) {
      std::cout << "<<< Test -1k >>> resolve boundary and step to layer "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // positive return: do the step to the layer
    step(state.stepping);
    state.stepping.stepSize = detail::ConstrainedStep(-beamPipeRadius);
    // step to the origin
    step(state.stepping);
    if (debug) {
      std::cout << "<<< Returned back at origin: "
                << toString(state.stepping.position()) << std::endl;
    }

    // ------ initialize testing --------------------
    if (debug) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<< INITIALIZE TESTING >>>>>>>>>>>>>>>>>>"
                << std::endl;
    }

    // re-initialize/update the stepping state
    state              = PropagatorState();
    state.stepping.pos = onBeamPipe;
    state.stepping.dir = momentum.unit();
    // initialize : should not return to stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // step size has been updated
    if (debug) {
      std::cout << "<<< Test 2a >>> initialize at BeamPipe  at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // update the stepping state
    state.stepping.pos = onBoundary;
    // initialize : should not return to stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // step size has been updated
    if (debug) {
      std::cout << "<<< Test 2b >>> initialize from Boundary at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // update the stepping state
    state.stepping.pos = on1stApproach;
    // initialize : should not return to stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // step size has been updated
    if (debug) {
      std::cout << "<<< Test 2c >>> initialize from Approach at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // update the stepping state
    state.stepping.pos = on1stSf1;
    // initialize : should not return to stepper
    BOOST_TEST(navigator.initialize(state, navCorr) == false);
    // step size has been updated
    if (debug) {
      std::cout << "<<< Test 2d >>> initialize from 1st Surface  at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }

    // ------ special case testing --------------------
    if (debug) {
      std::cout
          << "<<<<<<<<<<<<<<<<<<<<< SPECIAL CASE TESTING >>>>>>>>>>>>>>>>>>"
          << std::endl;
    }

    // Start from a boundary:
    // - try hitting layers, do not hit layers
    // - switch boundary to next and update step size to that
    // re-initialize/update the stepping state
    state               = PropagatorState();
    state.options.debug = debug;
    // let's shift the boundary position in z
    onBoundary[2] = 400.;
    momentum      = Vector3D(1., 1., 100.);
    // set the stepping parameters
    state.stepping.pos = onBoundary;
    state.stepping.dir = momentum.unit();
    // initialize : should not return to stepper
    BOOST_TEST(!navigator.initialize(state, navCorr) == true);
    // no surfaces to handle : do not return
    BOOST_TEST(navigator.handleSurfaces(state, navCorr) == false);
    // no layers to handle : do not return
    BOOST_TEST(navigator.handleLayers(state, navCorr) == false);
    // boundaries to handle : return
    BOOST_TEST(navigator.handleBoundaries(state, navCorr) == true);
    if (debug) {
      std::cout << "<<< Test 3a >>> start from boundary w/o layer hitting at "
                << toString(state.stepping.position()) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

}  // namespace Test

}  // namespace Acts
