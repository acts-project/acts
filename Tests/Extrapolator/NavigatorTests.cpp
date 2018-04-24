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

  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct Cache
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

    /// Navigation cache: the start surface
    const Surface* startSurface = nullptr;

    /// Navigation cache: the current surface
    const Surface* currentSurface = nullptr;

    /// Navigation cache: the target surface
    const Surface* targetSurface = nullptr;
    bool           targetReached = false;

    /// Debug output
    /// the string where things are stored (optionally)
    bool        debug        = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;
  };

  template <typename cache_t>
  NavigationParameters
  step(cache_t& cache)
  {
    // update the cache position
    cache.pos = cache.pos + cache.stepSize * cache.dir;
    // create navigation parameters
    return NavigationParameters(cache.pos, cache.dir);
  }

  // the surface cache & the creation of the geometry
  std::vector<std::unique_ptr<const Surface>> sCache;
  auto tGeometry = testGeometry<PlaneSurface>(sCache);

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
    Cache                  cache;
    cache.pos = position;
    cache.dir = momentum.unit();

    // (1) Initialization navigation from start point
    // - this will call resolveLayers() as well
    // - it will not return to the stepper
    BOOST_TEST(navigator.initialize(navPar, cache, result));
    // check that the currentVolume is set
    BOOST_TEST((result.currentVolume != nullptr));
    // check that the currentVolume is the startVolume
    BOOST_TEST((result.currentVolume == result.startVolume));
    // check that the currentSurface is rest to
    BOOST_TEST((cache.currentSurface == nullptr));
    // one layer has been found
    BOOST_TEST((result.navLayers.size() == 1));
    // the iterator should point to it
    BOOST_TEST((result.navLayerIter == result.navLayers.begin()));
    // cahce the beam pipe radius
    double beamPipeRadius = result.navLayerIter->intersection.position.perp();
    // step size has been updated
    BOOST_TEST((cache.stepSize == beamPipeRadius));
    if (navigator.debug) {
      std::cout << "<<< Test 1a >>> initialize at (0,0,0) " << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar             = step(cache);
    Vector3D onBeamPipe = navPar.position();

    // (2) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST(!(navigator.handeSurfaces(navPar, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handleLayers(navPar, cache, result)));
    BOOST_TEST((cache.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(navigator.handleBoundaries(navPar, cache, result));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1b >>> handleLayers, and set step to boundary  "
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar             = step(cache);
    Vector3D onBoundary = navPar.position();

    // (3) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST(!(navigator.handeSurfaces(navPar, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handleLayers(navPar, cache, result)));
    BOOST_TEST((cache.currentSurface != nullptr));
    // handleBoundaries should return true
    BOOST_TEST(navigator.handleBoundaries(navPar, cache, result));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1c >>> advance to boundary, initialize layers."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar                = step(cache);
    Vector3D on1stApproach = navPar.position();

    // (4) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST(!(navigator.handeSurfaces(navPar, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, cache, result)));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1d >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar           = step(cache);
    Vector3D on1stSf1 = navPar.position();

    // (5) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1e >>> advance to surface, update." << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (6) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1f >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (7) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1g >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (8) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1h >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (9) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST(!(navigator.handeSurfaces(navPar, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1i >>> advance to last surface, switch layer."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (10) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST(!(navigator.handeSurfaces(navPar, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handleLayers(navPar, cache, result)));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1j >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (5) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1k >>> advance to surface, update." << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // positive return: do the step
    navPar = step(cache);

    // (6) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // handeSurfaces should return false
    BOOST_TEST((navigator.handeSurfaces(navPar, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1l >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // Initialize navigation from beam pipe
    result = Navigator::result_type();
    // recreate cache
    cache     = Cache();
    cache.pos = onBeamPipe;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(cache.pos, momentum);
    BOOST_TEST(navPar.position() == onBeamPipe);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 2 >>> initialize at BeamPipe " << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // Initialize navigation from boundary surface
    result = Navigator::result_type();
    // recreate cache
    cache     = Cache();
    cache.pos = onBoundary;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(cache.pos, momentum);
    BOOST_TEST(navPar.position() == onBoundary);
    BOOST_TEST((navigator.initialize(navPar, cache, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 3 >>> initialize from Boundary " << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // Initialize navigation from approach surface of first layer
    result = Navigator::result_type();
    // recreate cache
    cache     = Cache();
    cache.pos = on1stApproach;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(cache.pos, momentum);
    BOOST_TEST(navPar.position() == on1stApproach);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 4 >>> initialize from Approach " << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }

    // Initialize navigation from 1st surface of
    result = Navigator::result_type();
    // recreate cache
    cache     = Cache();
    cache.pos = on1stSf1;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    navPar = NavigationParameters(cache.pos, momentum);
    BOOST_TEST(navPar.position() == on1stSf1);
    // this one avoids the stepping towards layer
    BOOST_TEST(!(navigator.initialize(navPar, cache, result)));
    // step size has been updated
    if (navigator.debug) {
      std::cout << "<<< Test 5 >>> initialize from 1st Surface " << std::endl;
      std::cout << cache.debugString << std::endl;
      cache.debugString = "";
    }
  }
}

}  // end of namespace Acts
