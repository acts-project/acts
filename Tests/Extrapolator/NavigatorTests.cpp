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
    NavigationDirection nav_dir = forward;

    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive sep size of the runge-kutta integration
    cstep step_size = 100 * units::_cm;

    /// Navigation cache: the start surface
    const Surface* start_surface = nullptr;

    /// Navigation cache: the current surface
    const Surface* current_surface = nullptr;

    /// Navigation cache: the target surface
    const Surface* target_surface = nullptr;
    bool           target_reached = false;

    /// Debug output
    /// the string where things are stored (optionally)
    bool        debug        = false;
    std::string debug_string = "";
    /// buffer & formatting for consistent output
    size_t debug_pfx_width = 30;
    size_t debug_msg_width = 50;
  };

  template <typename cache_t>
  NavigationParameters
  step(cache_t& cache)
  {
    // update the cache position
    cache.pos = cache.pos + cache.step_size * cache.dir;
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
    NavigationParameters nav_par(position, momentum);

    // the result and the cache
    Navigator::result_type result;
    Cache                  cache;
    cache.pos = position;
    cache.dir = momentum.unit();

    // (1) Initialization navigation from start point
    // - this will call resolve_layers() as well
    // - it will not return to the stepper
    BOOST_TEST(navigator.initialize(nav_par, cache, result));
    // check that the current_volume is set
    BOOST_TEST((result.current_volume != nullptr));
    // check that the current_volume is the start_volume
    BOOST_TEST((result.current_volume == result.start_volume));
    // check that the current_surface is rest to
    BOOST_TEST((cache.current_surface == nullptr));
    // one layer has been found
    BOOST_TEST((result.nav_layers.size() == 1));
    // the iterator should point to it
    BOOST_TEST((result.nav_layer_iter == result.nav_layers.begin()));
    // cahce the beam pipe radius
    double beamPipeRadius = result.nav_layer_iter->intersection.position.perp();
    // step size has been updated
    BOOST_TEST((cache.step_size == beamPipeRadius));
    if (navigator.debug) {
      std::cout << "<<< Test 1a >>> initialize at (0,0,0) " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par             = step(cache);
    Vector3D onBeamPipe = nav_par.position();

    // (2) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST(!(navigator.handle_surfaces(nav_par, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handle_layers(nav_par, cache, result)));
    BOOST_TEST((cache.current_surface != nullptr));
    // handle_boundaries should return true
    BOOST_TEST(navigator.handle_boundaries(nav_par, cache, result));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1b >>> handle_layers, and set step to boundary  "
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par             = step(cache);
    Vector3D onBoundary = nav_par.position();

    // (3) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST(!(navigator.handle_surfaces(nav_par, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST(!(navigator.handle_layers(nav_par, cache, result)));
    BOOST_TEST((cache.current_surface != nullptr));
    // handle_boundaries should return true
    BOOST_TEST(navigator.handle_boundaries(nav_par, cache, result));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1c >>> advance to boundary, initialize layers."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par                = step(cache);
    Vector3D on1stApproach = nav_par.position();

    // (4) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST(!(navigator.handle_surfaces(nav_par, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handle_layers(nav_par, cache, result)));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1d >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par           = step(cache);
    Vector3D on1stSf1 = nav_par.position();

    // (5) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST((navigator.handle_surfaces(nav_par, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1e >>> advance to surface, update." << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par           = step(cache);
    Vector3D on1stSf2 = nav_par.position();

    // (6) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST((navigator.handle_surfaces(nav_par, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1f >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par           = step(cache);
    Vector3D on1stSf3 = nav_par.position();

    // (7) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST((navigator.handle_surfaces(nav_par, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1g >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par           = step(cache);
    Vector3D on1stSf4 = nav_par.position();

    // (8) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST((navigator.handle_surfaces(nav_par, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1h >>> advance to next surface, update."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par           = step(cache);
    Vector3D on1stSf5 = nav_par.position();

    // (9) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST(!(navigator.handle_surfaces(nav_par, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handle_layers(nav_par, cache, result)));
    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1i >>> advance to last surface, switch layer."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // positive return: do the step
    nav_par                = step(cache);
    Vector3D on2ndApproach = nav_par.position();

    // (10) re-entering navigator:
    // initialize should return false
    BOOST_TEST(!(navigator.initialize(nav_par, cache, result)));
    // handle_surfaces should return false
    BOOST_TEST(!(navigator.handle_surfaces(nav_par, cache, result)));
    // handle_layer should sore the layer surface, but return false
    BOOST_TEST((navigator.handle_layers(nav_par, cache, result)));

    // the iterator should point to it
    if (navigator.debug) {
      std::cout << "<<< Test 1j >>> advance to layer, initialize surfaces."
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    /*
    // Initialize navigation from beam pipe
    // lets progess to the intersection with the beam-pipe
    Vector3D atBeampipe = result.nav_layer_iter->intersection.position;
    // the step size should be the radius of the beam pipe
    // clear result
    result = Navigator::result_type();
    // recreate cache
    cache     = Cache();
    cache.pos = atBeampipe;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    nav_par = NavigationParameters(cache.pos, momentum);
    BOOST_TEST(nav_par.position() != Vector3D(0., 0., 0.));

    // this will call resolve_layer - no returning to the stepper here
    BOOST_TEST(!navigator.initialize(nav_par, cache, result));
    BOOST_TEST(result.nav_layers.size() == 1);
    BOOST_TEST(cache.step_size == 0.);
    // the current surface should now be set
    BOOST_TEST((cache.current_surface != nullptr));

    if (navigator.debug) {
      std::cout << "<<< Test 3 >>> initialize at beam pipe " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // Now it's time to handle the current layers
    // we are already on the current layer
    BOOST_TEST(!navigator.handle_layers(nav_par, cache, result));
    // layers have to be reset
    BOOST_TEST(result.nav_layers.size() == 0);
    // The current surface should still point to the layer surface
    BOOST_TEST((cache.current_surface != nullptr));
    if (navigator.debug) {
      std::cout << "<<< Test 4 >>> handle layers " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // Let's test the exit of the BeamPipe volume
    BOOST_TEST(navigator.handle_boundaries(nav_par, cache, result));
    // We should have a boundary surface candidate
    BOOST_TEST(result.nav_boundaries.size() != 0.);
    // The step size should be updated towards the boundary
    BOOST_TEST(cache.step_size != 0.);
    if (navigator.debug) {
      std::cout << "<<< Test 5 >>> handle boundaries " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // Initialize navigation from boundary surface
    // advance to the boundary
    Vector3D atBoundary
        = nav_par.position() + cache.step_size * nav_par.momentum().unit();

    // recreate cache and nav pars
    cache     = Cache();
    cache.pos = atBoundary;
    cache.dir = momentum.unit();
    nav_par   = NavigationParameters(cache.pos, momentum);

    BOOST_TEST(navigator.handle_boundaries(nav_par, cache, result));
    if (navigator.debug) {
      std::cout << "<<< Test 6 >>> advance to boundary" << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // the step size should be the radius of the beam pipe
    // clear result
    result = Navigator::result_type();
    // recreate cache and nav pars
    cache     = Cache();
    cache.pos = atBoundary;
    cache.dir = momentum.unit();
    // recreate new navigation parameters
    BOOST_TEST(nav_par.position() == atBoundary);
    BOOST_TEST(navigator.initialize(nav_par, cache, result));

    if (navigator.debug) {
      std::cout << "<<< Test 7 >>> initialize at boundary " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // Now it's time to advance to first pixel layer
    Vector3D at1st = atBoundary + cache.step_size * nav_par.momentum().unit();
    cache.pos      = at1st;
    nav_par        = NavigationParameters(at1st, momentum);

    BOOST_TEST(navigator.handle_layers(nav_par, cache, result));

    if (navigator.debug) {
      std::cout << "<<< Test 8 >>> advance to and handle layers in pixels "
                << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // The starting point
    Vector3D atSf = at1st + cache.step_size * nav_par.momentum().unit();
    cache.pos = atSf;
    nav_par   = NavigationParameters(atSf, momentum);
    // and handle the surfaces
    bool progress = navigator.handle_surfaces(nav_par, cache, result);
    BOOST_TEST(progress);
    if (navigator.debug) {
      std::cout << "<<< Test 9 >>> advance to surface " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // clear result
    result = Navigator::result_type();
    // recreate cache and nav pars
    cache     = Cache();
    cache.pos = at1st;
    cache.dir = momentum.unit();
    nav_par   = NavigationParameters(at1st, momentum);

    // recreate new navigation parameters
    BOOST_TEST(nav_par.position() == at1st);
    BOOST_TEST(navigator.initialize(nav_par, cache, result));
    // we should have a current layer
    BOOST_TEST((result.start_layer!=nullptr));
    BOOST_TEST((result.nav_layer_iter->object!=result.start_layer));
    // the current layer should point to the begin
    if (navigator.debug) {
      std::cout << "<<< Test 10 >>> Initialize at 1st layer, approach surface "
    << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";

    }

    // clear result
    result = Navigator::result_type();
    // recreate cache and nav pars
    cache     = Cache();
    cache.pos = atSf;
    cache.dir = momentum.unit();
    nav_par   = NavigationParameters(atSf, momentum);

    // recreate new navigation parameters
    BOOST_TEST(nav_par.position() == atSf);
    BOOST_TEST(navigator.initialize(nav_par, cache, result));
    // we should have a current layer
    BOOST_TEST((result.start_layer!=nullptr));
    BOOST_TEST((result.nav_layer_iter->object!=result.start_layer));
    // the current layer should point to the begin
    if (navigator.debug) {
      std::cout << "<<< Test 11 >>> Initialize at 1st layer, sensor surface " <<
    std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";

    }

    */
  }
}

}  // end of namespace Acts
