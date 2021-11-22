// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Experimental/CylindricalContainerHelper.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/Tracer.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <vector>

#include "GeometryHelper.hpp"

namespace Acts {

using namespace UnitLiterals;

namespace Test {

// A minimal stapper for this test
struct Stepper {
  struct State {
    // Position
    Vector4 pos4 = Vector4(0., 0., 0., 0.);

    /// Direction
    Vector3 dir = Vector3(1., 0., 0.);
    /// Momentum
    ActsScalar p;

    /// Charge
    ActsScalar q;

    /// The navigation direction
    NavigationDirection navDir = forward;

    // Accummulated path length cache
    ActsScalar pathAccumulated = 0.;

    // Adaptive sep size of the runge-kutta integration
    ConstrainedStep stepSize = ConstrainedStep(100_cm);

    // A Geometry context
    GeometryContext geoContext = GeometryContext();
  };

  /// Global particle position accessor
  Vector3 position(const State& state) const {
    return state.pos4.segment<3>(Acts::ePos0);
  }

  /// Momentum direction accessor
  Vector3 direction(const State& state) const { return state.dir; }

  /// Momentum accessor
  ActsScalar momentum(const State& state) const { return state.p; }

  /// Charge access
  ActsScalar charge(const State& state) const { return state.q; }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
    detail::updateSingleStepSize<Stepper>(state, oIntersection, release);
  }
};

struct PropagatorState {
  Stepper::State stepping;
  Tracer::State navigation;
};

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(TracerStatusInitialize) {
  // Retrieve a detector instance
  auto detector = createDetector();

  Tracer::Config tConfig{};
  tConfig.world = detector.get();

  GeometryContext gctx = GeometryContext();

  // Central sector
  auto beamPipe = detector->lowest(gctx, Vector3(0., 0., 0.));
  auto layer0 = detector->lowest(gctx, Vector3(30., 0., 0.));
  auto gap = detector->lowest(gctx, Vector3(50., 0., 0.));
  auto layer1 = detector->lowest(gctx, Vector3(70., 0., 0.));

  // Create a simple stepper
  Stepper stepper;
  Tracer tracer(std::move(tConfig));

  // --------------------------------------------------------------------------
  PropagatorState pState;
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eUninitialized);

  // ------- First test ---------------------------
  // Status call from (0.,0.,0.) - initialize
  tracer.status(pState, stepper);
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsSurface);
  // No current surface set
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);

  // The beam pipe is the candidate
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  // ------- Second test ----------------------------
  // Status call from (23.,0.,0.) - initialize
  pState = PropagatorState();
  pState.stepping.pos4 = Vector4(23., 0., 0., 0.);
  tracer.status(pState, stepper);
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnSurface);
  // Current surface is set to the beampipe
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  // ------- Third test ----------------------------
  // Status call from (23.05,0.,0.) - initialize
  // A potential overstep tolerance will be ignored in an initialize call
  pState = PropagatorState();
  pState.stepping.pos4 = Vector4(23.05, 0., 0., 0.);
  tracer.status(pState, stepper);
  // There should be no surfaces
  BOOST_CHECK(pState.navigation.environment.surfaces.empty());
  // We should be on the way to a portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsPortal);
  // No current surface is set
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  // -------- Third test -----------------------------
  pState = PropagatorState();
  pState.stepping.pos4 = Vector4(27., 0., 0., 0.);

  // Status call from boundary (27.,0.,0.)
  tracer.status(pState, stepper);
  // The volume is layer0
  BOOST_CHECK(pState.navigation.environment.volume == layer0);
  // There should be surface candidates
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // There are 4 canidates (phi, a overlap at x,0,0)
  BOOST_CHECK(pState.navigation.environment.surfaces.size() == 4u);
  // Four portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 4u);
  // Current surface is an on portal surface
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);
  // It is the boundary surface of the beam pipe volume
  BOOST_CHECK(pState.navigation.currentSurface ==
              &(beamPipe->portals()[2]->surfaceRepresentation()));
  // We should be on the way to a portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnPortal);

  // Remember for the next step
  Vector3 atSurface =
      pState.navigation.environment.surfaces[0].intersection.position;

  // -------- Fourth test -----------------------------
  pState = PropagatorState();
  pState.stepping.pos4.block<3, 1>(0, 0) = atSurface;

  // Status call on a first surface
  tracer.status(pState, stepper);
  // The volume is layer0
  BOOST_CHECK(pState.navigation.environment.volume == layer0);
  // There should be surface candidates
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // There are 3 canidates (phi, a overlap at x,0,0) left
  BOOST_CHECK(pState.navigation.environment.surfaces.size() == 3u);
  // Four portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 4u);
  // Current surface is an on portal surface
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);
  // We should be on the way to a portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnSurface);

  // -------- Fifth test -----------------------------
  pState = PropagatorState();
  pState.stepping.pos4 = {
      80.,
      0.,
      0.,
      0.,
  };

  // Status call on a first surface
  tracer.status(pState, stepper);

  // End of world - no surfaces
  BOOST_CHECK(pState.navigation.environment.surfaces.empty());
  // End of world - no portals
  BOOST_CHECK(pState.navigation.environment.portals.empty());
  // End of world - no volume
  BOOST_CHECK(pState.navigation.environment.volume == nullptr);

  // --------------------------------------------------------------------------
  // Status calls from non-initialization
  pState = PropagatorState();
  tracer.status(pState, stepper);
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsSurface);

  // First test
  // Status call from (10.,0.,0.) - step status
  pState.stepping.pos4 = Vector4(10., 0., 0., 0.);
  tracer.status(pState, stepper);

  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsSurface);
  // Current surface is set to the beampipe
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  // Second test
  // Status call from beam pipe (23.,0.,0.) - step call
  pState.stepping.pos4 = Vector4(23., 0., 0., 0.);
  tracer.status(pState, stepper);

  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnSurface);
  // Current surface is set to the beampipe
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);

  // Third test
  // Status call from beam pipe (23.,0.,0.) - step call
  pState.stepping.pos4 = Vector4(23.05, 0., 0., 0.);
  tracer.status(pState, stepper);
  // During status call, there's no negative step size,
  // so the navigation is per se twoards the next portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsPortal);
  // Current surface is set to the beampipe
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);

  // Forth test
  // Status call from beam pipe volume boundary - step call
  pState.stepping.pos4 = Vector4(27., 0., 0., 0.);

  // Status call from boundary (27.,0.,0.)
  tracer.status(pState, stepper);
  // The volume is layer0
  BOOST_CHECK(pState.navigation.environment.volume == layer0);
  // There should be surface candidates
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // There are 4 canidates (phi, a overlap at x,0,0)
  BOOST_CHECK(pState.navigation.environment.surfaces.size() == 4u);
  // Four portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 4u);
  // Current surface is an on portal surface
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);
  // It is the boundary surface of the beam pipe volume
  BOOST_CHECK(pState.navigation.currentSurface ==
              &(beamPipe->portals()[2]->surfaceRepresentation()));
  // We should be on the way to a portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnPortal);

  // Fifth test
  pState.stepping.pos4.block<3, 1>(0, 0) = atSurface;

  // Status call on a first surface - step call
  tracer.status(pState, stepper);
  // The volume is layer0
  BOOST_CHECK(pState.navigation.environment.volume == layer0);
  // There should be surface candidates
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // There are 3 canidates (phi, a overlap at x,0,0), 
  // as one is excluded by on surface 
  BOOST_CHECK(pState.navigation.environment.surfaces.size() == 3u);
  // Four portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 4u);
  // Current surface is an on portal surface
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);
  // We should be on the way to a portal
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnSurface);
}

BOOST_AUTO_TEST_CASE(TracerNavigationChainNoOverstepping) {

  // Retrieve a detector instance
  auto detector = createDetector();

  Tracer::Config tConfig{};
  tConfig.world = detector.get();

  GeometryContext gctx = GeometryContext();

  // Central sector
  auto beamPipe = detector->lowest(gctx, Vector3(0., 0., 0.));
  auto layer0 = detector->lowest(gctx, Vector3(30., 0., 0.));
  auto gap = detector->lowest(gctx, Vector3(50., 0., 0.));
  auto layer1 = detector->lowest(gctx, Vector3(70., 0., 0.));

  // Create a simple stepper
  Stepper stepper;
  Tracer tracer(std::move(tConfig));

  // --------------------------------------------------------------------------
  PropagatorState pState;
  pState.stepping.dir = Vector3(1.,0.,0.);

  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eUninitialized);


  // ------- First test ---------------------------
  // Status call from (0.,0.,0.) - initialize
  tracer.status(pState, stepper);
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eTowardsSurface);
  // No current surface set
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);

  // The beam pipe is the candidate
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  // --- this would be the place for the actors to kick in

  // Target call 
  tracer.target(pState, stepper);
  // No current surface set
  BOOST_CHECK(pState.navigation.currentSurface == nullptr);
  // The beam pipe is the candidate
  BOOST_CHECK(not pState.navigation.environment.surfaces.empty());
  // Three portals
  BOOST_CHECK(pState.navigation.environment.portals.size() == 3u);
  // The volume is the beam pipe
  BOOST_CHECK(pState.navigation.environment.volume == beamPipe);

  CHECK_CLOSE_ABS(pState.stepping.stepSize, 23., s_onSurfaceTolerance);

  // Make the step towards the surface
  pState.stepping.pos4.block<3, 1>(0, 0) = pState.stepping.stepSize * pState.stepping.dir;

  tracer.status(pState, stepper);
  BOOST_CHECK(pState.navigation.environment.status ==
              DetectorEnvironment::eOnSurface);
  // Current surface is set to the beampipe
  BOOST_CHECK(pState.navigation.currentSurface != nullptr);


}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
