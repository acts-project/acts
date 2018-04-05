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
    cstep step_size = std::numeric_limits<double>::max();

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

  // the surface cache & the creation of the geometry
  std::vector<std::unique_ptr<const Surface>> sCache;
  auto                                        tGeometry = testGeometry(sCache);

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

    // Initialization navigation from start point
    BOOST_TEST(navigator.initialize(nav_par, cache, result));
    // check that the current_volume is set
    BOOST_TEST(result.current_volume != nullptr);
    BOOST_TEST(cache.current_surface == nullptr);
    BOOST_TEST(result.nav_layers.size() == 1);

    if (navigator.debug) {
      std::cout << "<<< Test 0 >>> initialize at (0,0,0) " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // this will be called by initialze() as well, doesn't
    // harm to call twice : initialze/resolve the layers
    BOOST_TEST(navigator.resolve_layers(nav_par, cache, result));

    if (navigator.debug) {
      std::cout << "<<< Test 1 >>> resolve_layers " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

    // Standard navigation:
    // we should have one layer now,
    // as there is only one material surface (it is the beam pipe)
    BOOST_TEST(result.nav_layers.size() == 1);
    // the nav_layer_iterator should point to the first element
    bool layer_set = (result.nav_layer_iter == result.nav_layers.begin());
    BOOST_TEST(layer_set);
    double beamPipeRadius = result.nav_layer_iter->intersection.position.perp();
    BOOST_TEST(cache.step_size == beamPipeRadius);

    // Also, if we try collect also the passive, we should only have one
    // layer -> navigation layers are always dark
    navigator.collectPassive = true;
    BOOST_TEST(navigator.resolve_layers(nav_par, cache, result));
    BOOST_TEST(result.nav_layers.size() == 1);
    // and the layer should still be set
    layer_set = (result.nav_layer_iter == result.nav_layers.begin());
    BOOST_TEST(layer_set);
    navigator.collectPassive = false;

    if (navigator.debug) {
      std::cout << "<<< Test 2 >>> resolve_layers with passive " << std::endl;
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }

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
    layer_set = (result.nav_layer_iter == result.nav_layers.begin());
    BOOST_TEST(layer_set);
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
    Vector3D atSf = at1st;

    if (navigator.debug) {
      std::cout << "<<< Test 9 >>> advance to and handle surface in pixels "
                << std::endl;
    }

    // bool progress = true;
    // while (progress && result.nav_surface_iter != result.nav_surfaces.end()){

    // progress to the surface
    atSf += cache.step_size * nav_par.momentum().unit();
    cache.pos = atSf;
    nav_par   = NavigationParameters(atSf, momentum);
    // and handle the surfaces
    bool progress = navigator.handle_surfaces(nav_par, cache, result);
    BOOST_TEST(progress);
    if (navigator.debug) {
      std::cout << cache.debug_string << std::endl;
      cache.debug_string = "";
    }
  }
}

}  // end of namespace Acts
