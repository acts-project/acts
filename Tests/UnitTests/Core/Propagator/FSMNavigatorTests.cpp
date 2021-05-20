// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/FSMNavigator.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

using namespace Acts::UnitLiterals;

namespace Acts {

using VectorHelpers::perp;

namespace Test {
GeometryContext tgContext = GeometryContext();

/// This is a simple cache struct to mimic a Stepper
struct PseudoStepper {
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
    double p;

    /// Charge
    double q;

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

  void step(State& sstate, double fraction) {
    // update the cache position
    double ssize = sstate.stepSize * fraction;
    std::cout << "PseudoStepper: Performing step with size: " << ssize
              << " along [" << sstate.dir.transpose() << "]: " << std::endl;
    Vector4 prev = sstate.pos4;
    sstate.pos4.head<3>() += (sstate.stepSize * fraction) * sstate.dir;
    auto rz = [](const Vector4& v) -> std::string {
      return std::to_string(VectorHelpers::perp(v)) + "," +
             std::to_string(v[eFreePos2]);
    };
    std::cout << "               [" << prev.transpose();
    std::cout << "] -> [" << sstate.pos4.transpose() << "]" << std::endl;

    std::cout << "               [" << rz(prev);
    std::cout << "] -> [" << rz(sstate.pos4.transpose()) << "]" << std::endl;

    // create navigation parameters
    return;
  }

  /// State resetter
  void resetState(State& /*unused*/, const BoundVector& /*unused*/,
                  const BoundSymMatrix& /*unused*/, const Surface& /*unused*/,
                  const NavigationDirection /*unused*/,
                  const double /*unused*/) const {}

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

  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck) const {
    return detail::updateSingleSurfaceStatus<PseudoStepper>(*this, state,
                                                            surface, bcheck);
  }

  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
    detail::updateSingleStepSize<PseudoStepper>(state, oIntersection, release);
  }

  void setStepSize(State& state, double stepSize,
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

  Result<BoundState> boundState(State& state, const Surface& surface,
                                bool /*unused*/
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

  CurvilinearState curvilinearState(State& state, bool /*unused*/
  ) const {
    CurvilinearTrackParameters parameters(state.pos4, state.dir, state.p,
                                          state.q);
    // Create the bound state
    CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                               state.pathAccumulated};
    return curvState;
  }

  void update(State& /*state*/, const FreeVector& /*pars*/,
              const Covariance& /*cov*/) const {}

  void update(State& /*state*/, const Vector3& /*uposition*/,
              const Vector3& /*udirection*/, double /*up*/,
              double /*time*/) const {}

  void transportCovarianceToCurvilinear(State& /*state*/) const {}

  void transportCovarianceToBound(State& /*unused*/,
                                  const Surface& /*surface*/) const {}

  Vector3 getField(State& /*state*/, const Vector3& /*pos*/) const {
    // get the field from the cell
    return Vector3(0., 0., 0.);
  }
};

static_assert(StepperConcept<PseudoStepper>,
              "Dummy stepper does not fulfill concept");

/// This is a simple cache struct to mimic the
/// Propagator cache
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

    LoggerWrapper logger{getDummyLogger()};
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

  /// The Stepper state - internal state of the Stepper
  PseudoStepper::State stepping;

  /// Navigation state - internal state of the Navigator
  FSMNavigator::State navigation;

  // The context cache for this propagation
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

  Vector4 position(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  PseudoStepper stepper;
  PropagatorState state;

  auto step = [&](double fraction = 1.) {
    stepper.step(state.stepping, fraction);
    navigator.status(state, stepper);
    std::cout << "ACTORS / ABORTERS" << std::endl;
    navigator.target(state, stepper);
  };

  state.stepping.pos4 = position;
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
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  step(1.0);
  // we're on the beampipe layer now
  BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                    GeometryIdentifier{}.setVolume(2).setLayer(2));

  step(0.5);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  step(1.0);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                    GeometryIdentifier{}.setVolume(3).setBoundary(4));

  step(0.5);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  step(1.0);
  BOOST_CHECK_EQUAL(
      state.navigation.currentSurface->geometryId(),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setApproach(1));

  step(0.5);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  step(1.0);

  // // we're on the first surface now
  BOOST_CHECK_EQUAL(
      state.navigation.currentSurface->geometryId(),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(122));

  step(0.5);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  step(1.0);

  // we're on the second surface now
  BOOST_CHECK_EQUAL(
      state.navigation.currentSurface->geometryId(),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(123));
}

// the debug boolean
bool debug = true;

/// @brief Method for testing vectors in @c Navigator::State
///
/// @param [in] state Navigator state
/// @param [in] navSurf Number of navigation surfaces
/// @param [in] navLay Number of navigation layers
/// @param [in] navBound Number of navigation boundaries
/// @param [in] extSurf Number of external surfaces
bool testNavigatorStateVectors(FSMNavigator::State& state, size_t navSurf,
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
    FSMNavigator::State& state, const TrackingVolume* worldVol,
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

// BOOST_AUTO_TEST_CASE(FSMNavigator_status_methods) {
//   // create a navigator
//   FSMNavigator::Config cfg;
//   cfg.resolveSensitive = false;
//   cfg.resolveMaterial = false;
//   cfg.resolvePassive = false;
//   FSMNavigator navigator{cfg};

//   // position and direction vector
//   Vector4 position4(0., 0., 0, 0);
//   Vector3 momentum(1., 1., 0);

//   // the propagator cache
//   PropagatorState state;
//   state.options.debug = debug;

//   // the stepper cache
//   state.stepping.pos4 = position4;
//   state.stepping.dir = momentum.normalized();

//   // Stepper
//   PseudoStepper stepper;

//   //
//   // (1) Test for inactivity
//   //
//   // Run without anything present
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr));

//   // Run with geometry but without resolving
//   cfg.trackingGeometry = tGeometry;
//   navigator = FSMNavigator{cfg};
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr));

//   // Run with geometry and resolving but broken navigation for various
//   reasons cfg.resolveSensitive = true; cfg.resolveMaterial = true;
//   cfg.resolvePassive = true;
//   navigator = FSMNavigator{cfg};
//   state.navigation.navigationBreak = true;
//   // a) Because target is reached
//   state.navigation.targetReached = true;
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr));
//   // b) Beacause of no target surface
//   state.navigation.targetReached = false;
//   state.navigation.targetSurface = nullptr;
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr, nullptr,
//                                          nullptr, nullptr, nullptr));
//   // c) Because the target surface is reached
//   const Surface* startSurf = tGeometry->getBeamline();
//   state.stepping.pos4.segment<3>(Acts::ePos0) =
//       startSurf->center(state.geoContext);
//   const Surface* targetSurf = startSurf;
//   state.navigation.targetSurface = targetSurf;
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, nullptr, nullptr,
//                                          nullptr, nullptr, targetSurf,
//                                          nullptr, nullptr, nullptr,
//                                          targetSurf));

//   //
//   // (2) Test the initialisation
//   //
//   // a) Initialise without additional information
//   state.navigation = FSMNavigator::State();
//   state.stepping.pos4 << 0., 0., 0., 0.;
//   const TrackingVolume* worldVol = tGeometry->highestTrackingVolume();
//   const TrackingVolume* startVol = tGeometry->lowestTrackingVolume(
//       state.geoContext, stepper.position(state.stepping));
//   const Layer* startLay = startVol->associatedLayer(
//       state.geoContext, stepper.position(state.stepping));
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, worldVol,
//   startVol,
//                                          startLay, nullptr, nullptr,
//                                          startVol, nullptr, nullptr,
//                                          nullptr));

//   // b) Initialise having a start surface
//   state.navigation = FSMNavigator::State();
//   state.navigation.startSurface = startSurf;
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, worldVol,
//   startVol,
//                                          startLay, startSurf, startSurf,
//                                          startVol, nullptr, nullptr,
//                                          nullptr));

//   // c) Initialise having a start volume
//   state.navigation = FSMNavigator::State();
//   state.navigation.startVolume = startVol;
//   navigator.status(state, stepper);
//   BOOST_CHECK(testNavigatorStateVectors(state.navigation, 0u, 0u, 0u, 0u));
//   BOOST_CHECK(testNavigatorStatePointers(state.navigation, worldVol,
//   startVol,
//                                          startLay, nullptr, nullptr,
//                                          startVol, nullptr, nullptr,
//                                          nullptr));
// }

}  // namespace Test
}  // namespace Acts
