// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"

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
    using BoundState = std::tuple<BoundParameters, Jacobian, double>;
    using CurvilinearState =
        std::tuple<CurvilinearParameters, Jacobian, double>;
    using BField = int;

    template <typename, typename>
    using return_parameter_type = void;

    /// This is a simple cache struct to mimic the
    /// Stepper cache in the propagation
    struct State {
      /// Position
      Vector3D pos = Vector3D(0., 0., 0.);

      /// Direction
      Vector3D dir = Vector3D(1., 0., 0.);

      /// Momentum
      double p;

      /// Charge
      double q;

      /// Time
      double t;

      /// the navigation direction
      NavigationDirection navDir = forward;

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

    /// Global particle position accessor
    Vector3D position(const State& state) const { return state.pos; }

    /// Momentum direction accessor
    Vector3D direction(const State& state) const { return state.dir; }

    /// Momentum accessor
    double momentum(const State& state) const { return state.p; }

    /// Charge access
    double charge(const State& state) const { return state.q; }

    /// Time access
    double time(const State& state) const { return state.t; }

    /// Overstep limit access
    double overstepLimit(const State& /*state*/) const {
      return s_onSurfaceTolerance;
    }

    Intersection::Status updateSurfaceStatus(
        State& state, const Surface& surface,
        const BoundaryCheck& bcheck) const {
      return detail::updateSingleSurfaceStatus<Stepper>(*this, state, surface,
                                                        bcheck);
    }

    template <typename object_intersection_t>
    void updateStepSize(State& state,
                        const object_intersection_t& oIntersection,
                        bool release = true) const {
      detail::updateSingleStepSize<Stepper>(state, oIntersection, release);
    }

    void setStepSize(
        State& state, double stepSize,
        ConstrainedStep::Type stype = ConstrainedStep::actor) const {
      state.previousStepSize = state.stepSize;
      state.stepSize.update(stepSize, stype, true);
    }

    void releaseStepSize(State& state) const {
      state.stepSize.release(ConstrainedStep::actor);
    }

    std::string outputStepSize(const State& state) const {
      return state.stepSize.toString();
    }

    BoundState boundState(State& state, const Surface& surface) const {
      BoundParameters parameters(tgContext, std::nullopt, state.pos,
                                 state.p * state.dir, state.q, state.t,
                                 surface.getSharedPtr());
      BoundState bState{std::move(parameters), Jacobian::Identity(),
                        state.pathAccumulated};
      return bState;
    }

    CurvilinearState curvilinearState(State& state) const {
      CurvilinearParameters parameters(std::nullopt, state.pos,
                                       state.p * state.dir, state.q, state.t);
      // Create the bound state
      CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                                 state.pathAccumulated};
      return curvState;
    }

    void update(State& /*state*/, const BoundParameters& /*pars*/) const {}

    void update(State& /*state*/, const Vector3D& /*uposition*/,
                const Vector3D& /*udirection*/, double /*up*/,
                double /*time*/) const {}

    void covarianceTransport(State& /*state*/) const {}

    void covarianceTransport(State& /*unused*/,
                             const Surface& /*surface*/) const {}

    Vector3D getField(State& /*state*/, const Vector3D& /*pos*/) const {
      // get the field from the cell
      return Vector3D(0., 0., 0.);
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
  sstate.pos = sstate.pos + sstate.stepSize * sstate.dir;
  // create navigation parameters
  return;
}

// the surface cache & the creation of the geometry

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// the debug boolean
bool debug = true;

BOOST_AUTO_TEST_CASE(Navigator_methods) {
  // create a navigator
  Navigator navigator;
  navigator.trackingGeometry = tGeometry;
  navigator.resolveSensitive = true;
  navigator.resolveMaterial = true;
  navigator.resolvePassive = false;

  // position and direction vector
  Vector3D position(0., 0., 0);
  Vector3D momentum(1., 1., 0);

  // the propagator cache
  PropagatorState state;
  state.options.debug = debug;

  // the stepper cache
  state.stepping.pos = position;
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
  // The iterator should points to the begin
  BOOST_CHECK(state.navigation.navLayerIter ==
              state.navigation.navLayers.begin());
  // Cache the beam pipe radius
  double beamPipeR = perp(state.navigation.navLayerIter->intersection.position);
  // step size has been updated
  CHECK_CLOSE_ABS(state.stepping.stepSize, beamPipeR, s_onSurfaceTolerance);
  if (debug) {
    std::cout << "<<< Test 1a >>> initialize at "
              << toString(state.stepping.pos) << std::endl;
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
  // The iterator still points to the begin
  BOOST_CHECK(
      (state.navigation.navLayerIter == state.navigation.navLayers.begin()));
  // ACTORS-ABORTERS-TARGET
  navigator.target(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1b >>> step to the BeamPipe at  "
              << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
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
                << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
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
                << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
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
                << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
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
                << toString(state.stepping.pos) << std::endl;
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
              << toString(state.stepping.pos) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }
}

}  // namespace Test
}  // namespace Acts
