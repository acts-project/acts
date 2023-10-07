// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CubicBVHTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <memory>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::perp;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

/// This is a simple cache struct to mimic the
/// Propagator cache
struct PropagatorState {
  /// This is a simple cache struct to mimic a Stepper
  struct Stepper {
    // comply with concept
    using Jacobian = BoundMatrix;
    using Covariance = BoundSymMatrix;
    using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
    using CurvilinearState =
        std::tuple<CurvilinearTrackParameters, Jacobian, double>;
    using BField = int;

    template <typename, typename>
    using return_parameter_type = void;

    /// This is a simple cache struct to mimic the
    /// Stepper cache in the propagation
    struct State {
      /// Position
      Vector4 pos4 = Vector4(0., 0., 0., 0.);

      /// Direction
      Vector3 dir = Vector3(1., 0., 0.);

      /// Momentum
      double p = 0;

      /// Charge
      double q = 0;

      /// the navigation direction
      NavigationDirection navDir = NavigationDirection::Forward;

      // accummulated path length cache
      double pathAccumulated = 0.;

      // adaptive sep size of the runge-kutta integration
      ConstrainedStep stepSize = ConstrainedStep(100_cm);

      // Previous step size for overstep estimation (ignored here)
      double previousStepSize = 0.;

      /// The tolerance for the stepping
      double tolerance = s_onSurfaceTolerance;

      GeometryContext geoContext = GeometryContext();
    };

    /// State resetter
    void resetState(State& /*state*/, const BoundVector& /*boundParams*/,
                    const BoundSymMatrix& /*cov*/, const Surface& /*surface*/,
                    const NavigationDirection /*navDir*/,
                    const double /*stepSize*/) const {}

    /// Global particle position accessor
    Vector3 position(const State& state) const {
      return state.pos4.segment<3>(Acts::ePos0);
    }

    /// Time access
    double time(const State& state) const { return state.pos4[Acts::eTime]; }

    /// Momentum direction accessor
    Vector3 direction(const State& state) const { return state.dir; }

    /// Momentum accessor
    double momentum(const State& state) const { return state.p; }

    /// Charge access
    double charge(const State& state) const { return state.q; }

    /// Overstep limit access
    double overstepLimit(const State& /*state*/) const {
      return s_onSurfaceTolerance;
    }

    Intersection3D::Status updateSurfaceStatus(State& state,
                                               const Surface& surface,
                                               const BoundaryCheck& bcheck,
                                               const Logger& logger) const {
      return detail::updateSingleSurfaceStatus<Stepper>(*this, state, surface,
                                                        bcheck, logger);
    }

    template <typename object_intersection_t>
    void updateStepSize(State& state,
                        const object_intersection_t& oIntersection,
                        bool release = true) const {
      detail::updateSingleStepSize<Stepper>(state, oIntersection, release);
    }

    void setStepSize(State& state, double stepSize,
                     ConstrainedStep::Type stype = ConstrainedStep::actor,
                     bool release = true) const {
      state.previousStepSize = state.stepSize.value();
      state.stepSize.update(stepSize, stype, release);
    }

    double getStepSize(const State& state, ConstrainedStep::Type stype) const {
      return state.stepSize.value(stype);
    }

    void releaseStepSize(State& state) const {
      state.stepSize.release(ConstrainedStep::actor);
    }

    std::string outputStepSize(const State& state) const {
      return state.stepSize.toString();
    }

    Result<BoundState> boundState(
        State& state, const Surface& surface, bool /*transportCov*/,
        const FreeToBoundCorrection& /*freeToBoundCorrection*/
    ) const {
      auto bound =
          BoundTrackParameters::create(surface.getSharedPtr(), tgContext,
                                       state.pos4, state.dir, state.p, state.q);
      if (!bound.ok()) {
        return bound.error();
      }
      BoundState bState{std::move(*bound), Jacobian::Identity(),
                        state.pathAccumulated};
      return bState;
    }

    CurvilinearState curvilinearState(State& state, bool /*transportCov*/
    ) const {
      CurvilinearTrackParameters parameters(state.pos4, state.dir, state.p,
                                            state.q);
      // Create the bound state
      CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                                 state.pathAccumulated};
      return curvState;
    }

    void update(State& /*state*/, const FreeVector& /*freePars*/,
                const BoundVector& /*boundPars*/, const Covariance& /*cov*/,
                const Surface& /*surface*/) const {}

    void update(State& /*state*/, const Vector3& /*uposition*/,
                const Vector3& /*udirection*/, double /*up*/,
                double /*time*/) const {}

    void transportCovarianceToCurvilinear(State& /*state*/) const {}

    void transportCovarianceToBound(
        State& /*state*/, const Surface& /*surface*/,
        const FreeToBoundCorrection& /*freeToBoundCorrection*/) const {}

    Result<Vector3> getField(State& /*state*/, const Vector3& /*pos*/) const {
      // get the field from the cell
      return Result<Vector3>::success({0., 0., 0.});
    }
  };

  static_assert(StepperConcept<Stepper>,
                "Dummy stepper does not fulfill concept");

  /// emulate the options template
  struct Options {
    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;

    const Acts::Logger& logger = Acts::getDummyLogger();
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
  Stepper::State stepping;

  /// Navigation state - internal state of the Navigator
  Navigator::State navigation;

  // The context cache for this propagation
  GeometryContext geoContext = GeometryContext();
};

template <typename stepper_state_t>
void step(stepper_state_t& sstate) {
  // update the cache position
  sstate.pos4[Acts::ePos0] += sstate.stepSize.value() * sstate.dir[Acts::eMom0];
  sstate.pos4[Acts::ePos1] += sstate.stepSize.value() * sstate.dir[Acts::eMom1];
  sstate.pos4[Acts::ePos2] += sstate.stepSize.value() * sstate.dir[Acts::eMom2];
  // create navigation parameters
  return;
}

/// @brief Method for testing vectors in @c Navigator::State
///
/// @param [in] state Navigator state
/// @param [in] navSurf Number of navigation surfaces
/// @param [in] navLay Number of navigation layers
/// @param [in] navBound Number of navigation boundaries
/// @param [in] extSurf Number of external surfaces
bool testNavigatorStateVectors(Navigator::State& state, size_t navSurf,
                               size_t navLay, size_t navBound, size_t extSurf) {
  return ((state.navSurfaces.size() == navSurf) &&
          (state.navLayers.size() == navLay) &&
          (state.navBoundaries.size() == navBound) &&
          (state.externalSurfaces.size() == extSurf));
}

/// @brief Method for testing pointers in @c Navigator::State
///
/// @param [in] state Navigation state
/// @param [in] worldVol World volume
/// @param [in] startVol Start volume
/// @param [in] startLay Start layer
/// @param [in] startSurf Start surface
/// @param [in] currSurf Current surface
/// @param [in] currVol Current volume
/// @param [in] targetVol Target volume
/// @param [in] targetLay Target layer
/// @param [in] targetSurf Target surface
bool testNavigatorStatePointers(
    Navigator::State& state, const TrackingVolume* worldVol,
    const TrackingVolume* startVol, const Layer* startLay,
    const Surface* startSurf, const Surface* currSurf,
    const TrackingVolume* currVol, const TrackingVolume* targetVol,
    const Layer* targetLay, const Surface* targetSurf) {
  return (
      (state.worldVolume == worldVol) && (state.startVolume == startVol) &&
      (state.startLayer == startLay) && (state.startSurface == startSurf) &&
      (state.currentSurface == currSurf) && (state.currentVolume == currVol) &&
      (state.targetVolume == targetVol) && (state.targetLayer == targetLay) &&
      (state.targetSurface == targetSurf));
}
// the surface cache & the creation of the geometry

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// the debug boolean
bool debug = true;

BOOST_AUTO_TEST_CASE(Navigator_status_methods) {
  // position and direction vector
  Vector4 position4(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  // the propagator cache
  PropagatorState state;
  state.options.debug = debug;

  // the stepper cache
  state.stepping.pos4 = position4;
  state.stepping.dir = momentum.normalized();

  // Stepper
  PropagatorState::Stepper stepper;

  //
  // (1) Test for inactivity
  //
  // Run without anything present
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    Navigator navigator{navCfg};

    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
                                           nullptr, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));
  }

  // Run with geometry but without resolving
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
                                           nullptr, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));
  }

  // Run with geometry and resolving but broken navigation for various reasons
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = true;
    navCfg.resolveMaterial = true;
    navCfg.resolvePassive = true;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    state.navigation.navigationBreak = true;
    // a) Because target is reached
    state.navigation.targetReached = true;
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
                                           nullptr, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));

    // b) Beacause of no target surface
    state.navigation.targetReached = false;
    state.navigation.targetSurface = nullptr;
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
                                           nullptr, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));
    // c) Because the target surface is reached
    const Surface* startSurf = tGeometry->getBeamline();
    state.stepping.pos4.segment<3>(Acts::ePos0) =
        startSurf->center(state.geoContext);
    const Surface* targetSurf = startSurf;
    state.navigation.targetSurface = targetSurf;
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(
        state.navigation, nullptr, nullptr, nullptr, nullptr, targetSurf,
        nullptr, nullptr, nullptr, targetSurf));

    //
    // (2) Test the initialisation
    //
    // a) Initialise without additional information
    state.navigation = Navigator::State();
    state.stepping.pos4 << 0., 0., 0., 0.;
    const TrackingVolume* worldVol = tGeometry->highestTrackingVolume();
    const TrackingVolume* startVol = tGeometry->lowestTrackingVolume(
        state.geoContext, stepper.position(state.stepping));
    const Layer* startLay = startVol->associatedLayer(
        state.geoContext, stepper.position(state.stepping));
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, worldVol, startVol,
                                           startLay, nullptr, nullptr, startVol,
                                           nullptr, nullptr, nullptr));

    // b) Initialise having a start surface
    state.navigation = Navigator::State();
    state.navigation.startSurface = startSurf;
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(
        state.navigation, worldVol, startVol, startLay, startSurf, startSurf,
        startVol, nullptr, nullptr, nullptr));

    // c) Initialise having a start volume
    state.navigation = Navigator::State();
    state.navigation.startVolume = startVol;
    navigator.status(state, stepper);
    BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state.navigation, worldVol, startVol,
                                           startLay, nullptr, nullptr, startVol,
                                           nullptr, nullptr, nullptr));
  }
}

BOOST_AUTO_TEST_CASE(Navigator_target_methods) {
  // create a navigator
  Navigator::Config navCfg;
  navCfg.trackingGeometry = tGeometry;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator{navCfg};

  // create a navigator for the Bounding Volume Hierarchy test
  CubicBVHTrackingGeometry grid(20, 1000, 5);
  Navigator::Config bvhNavCfg;
  bvhNavCfg.trackingGeometry = grid.trackingGeometry;
  bvhNavCfg.resolveSensitive = true;
  bvhNavCfg.resolveMaterial = true;
  bvhNavCfg.resolvePassive = false;
  Navigator BVHNavigator{bvhNavCfg};

  // position and direction vector
  Vector4 position4(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  // the propagator cache
  PropagatorState state;
  state.options.debug = debug;

  // the stepper cache
  state.stepping.pos4 = position4;
  state.stepping.dir = momentum.normalized();

  // foward navigation ----------------------------------------------
  if (debug) {
    std::cout << "<<<<<<<<<<<<<<<<<<<<< FORWARD NAVIGATION >>>>>>>>>>>>>>>>>>"
              << std::endl;
  }

  // Stepper
  PropagatorState::Stepper stepper;

  // (1) Initialization navigation from start point
  // - this will call resolveLayers() as well
  // - and thus should call a return to the stepper
  navigator.status(state, stepper);
  // Check that the currentVolume is set
  BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
  // Check that the currentVolume is the startVolume
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  // Check that the currentSurface is reset to:
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  // No layer has been found
  BOOST_CHECK_EQUAL(state.navigation.navLayers.size(), 0u);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);
  // A layer has been found
  BOOST_CHECK_EQUAL(state.navigation.navLayers.size(), 1u);
  // The index should points to the begin
  BOOST_CHECK(state.navigation.navLayerIndex == 0);
  // Cache the beam pipe radius
  double beamPipeR = perp(state.navigation.navLayer().intersection.position);
  // step size has been updated
  CHECK_CLOSE_ABS(state.stepping.stepSize.value(), beamPipeR,
                  s_onSurfaceTolerance);
  if (debug) {
    std::cout << "<<< Test 1a >>> initialize at "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    // Clear the debug string for the next test
    state.options.debugString = "";
  }

  // Do the step towards the beam pipe
  step(state.stepping);

  // (2) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // Check that the currentVolume is the still startVolume
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  // The layer number has not changed
  BOOST_CHECK_EQUAL(state.navigation.navLayers.size(), 1u);
  // The index still points to the begin
  BOOST_CHECK(state.navigation.navLayerIndex == 0);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1b >>> step to the BeamPipe at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Do the step towards the boundary
  step(state.stepping);

  // (3) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1c >>> step to the Boundary at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // positive return: do the step
  step(state.stepping);
  // (4) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1d >>> step to 1st layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on first layer
  for (size_t isf = 0; isf < 5; ++isf) {
    step(state.stepping);
    // (5-9) re-entering navigator:
    // STATUS
    navigator.status(state, stepper);
    // ACTORS-ABORTERS-TARGET
    navigator.target(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1e-1i >>> step within 1st layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (10) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1j >>> step to 2nd layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on second layer
  for (size_t isf = 0; isf < 5; ++isf) {
    step(state.stepping);
    // (11-15) re-entering navigator:
    // STATUS
    navigator.status(state, stepper);
    // ACTORS-ABORTERS-TARGET
    navigator.target(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1k-1o >>> step within 2nd layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (16) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1p >>> step to 3rd layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on third layer
  for (size_t isf = 0; isf < 3; ++isf) {
    step(state.stepping);
    // (17-19) re-entering navigator:
    // STATUS
    navigator.status(state, stepper);
    // ACTORS-ABORTERS-TARGET
    navigator.target(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1q-1s >>> step within 3rd layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (20) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1t >>> step to 4th layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on second layer
  for (size_t isf = 0; isf < 3; ++isf) {
    step(state.stepping);
    // (21-23) re-entering navigator:
    // STATUS
    navigator.status(state, stepper);
    // ACTORS-ABORTERS-TARGET
    navigator.target(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1t-1v >>> step within 4th layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (24) re-entering navigator:
  // STATUS
  navigator.status(state, stepper);
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1w >>> step to boundary at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // test the navigation in a bounding volume hierarchy
  // ----------------------------------------------
  if (debug) {
    std::cout << "<<<<<<<<<<<<<<<<<<<<< BVH Navigation >>>>>>>>>>>>>>>>>>"
              << std::endl;
  }

  // position and direction vector
  Vector4 BVHPosition4(0., 0., 0, 0);
  Vector3 BVHMomentum(1., 1., 0.);

  // the propagator cache
  PropagatorState BVHState;
  BVHState.options.debug = debug;

  // the stepper cache
  BVHState.stepping.pos4 = BVHPosition4;
  BVHState.stepping.dir = BVHMomentum.normalized();

  // Stepper
  PropagatorState::Stepper BVHStepper;

  BVHNavigator.status(BVHState, BVHStepper);

  // Check that the currentVolume is set
  BOOST_CHECK_NE(BVHState.navigation.currentVolume, nullptr);
  // Check that the currentVolume is the startVolume
  BOOST_CHECK_EQUAL(BVHState.navigation.currentVolume,
                    BVHState.navigation.startVolume);
  // Check that the currentSurface is reset to:
  BOOST_CHECK_EQUAL(BVHState.navigation.currentSurface, nullptr);
  // No layer has been found
  BOOST_CHECK_EQUAL(BVHState.navigation.navLayers.size(), 0u);
  // ACTORS-ABORTERS-TARGET
  navigator.target(BVHState, BVHStepper);
  // Still no layer
  BOOST_CHECK_EQUAL(BVHState.navigation.navLayers.size(), 0u);
  // Surfaces have been found
  BOOST_CHECK_EQUAL(BVHState.navigation.navSurfaces.size(), 42u);
}

}  // namespace Test
}  // namespace Acts
