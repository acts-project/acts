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
#include "Acts/Propagator/Navigator.hpp"
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
  Navigator::State navigation;

  // The context cache for this propagation
  GeometryContext geoContext = GeometryContext();
};

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

BOOST_AUTO_TEST_CASE(Navigation) {
  // FSMNavigator::Config cfg;
  // cfg.trackingGeometry = tGeometry;
  // cfg.resolveSensitive = true;
  // cfg.resolveMaterial = true;
  // cfg.resolvePassive = true;
  Navigator navigator(tGeometry);

  Vector4 position(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  CurvilinearTrackParameters params{position, momentum.normalized(),
                                    momentum.norm(), 1};

  const auto* startSurface = &params.referenceSurface();

  PseudoStepper stepper;
  PropagatorState state;

  state.navigation.startSurface = startSurface;

  auto logger = getDefaultLogger("Navigator", Logging::Level::VERBOSE);
  state.options.logger = LoggerWrapper{*logger};

  auto step = [&](double fraction = 1.) {
    stepper.step(state.stepping, fraction);
    // navigator.status(state, stepper);
    // std::cout << "ACTORS / ABORTERS" << std::endl;
    // navigator.target(state, stepper);
  };

  auto status = [&]() {
    std::cout << "STATUS" << std::endl;
    navigator.status(state, stepper);
  };

  auto target = [&]() {
    std::cout << "TARGET" << std::endl;
    navigator.target(state, stepper);
  };

  state.stepping.pos4 = position;
  state.stepping.dir = momentum.normalized();

  status();

  // currentVolume has been set
  BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, startSurface);
  BOOST_CHECK_EQUAL(state.navigation.startSurface, startSurface);

  BOOST_CHECK_EQUAL(state.navigation.currentVolume->volumeName(),
                    "BeamPipe::Barrel");

  // status is done, we should still be in initial state
  // propagator now calls target
  target();

  std::vector<GeometryIdentifier> surfaceSequence{
      GeometryIdentifier{}.setVolume(2).setLayer(2),
      GeometryIdentifier{}.setVolume(3).setBoundary(4),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(122),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(123),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(106),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setSensitive(107),
      GeometryIdentifier{}.setVolume(3).setLayer(2).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(244),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(212),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(245),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setSensitive(213),
      GeometryIdentifier{}.setVolume(3).setLayer(4).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setSensitive(397),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setSensitive(345),
      GeometryIdentifier{}.setVolume(3).setLayer(6).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setApproach(1),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setSensitive(595),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setSensitive(517),
      GeometryIdentifier{}.setVolume(3).setLayer(8).setApproach(2),
      GeometryIdentifier{}.setVolume(3).setBoundary(3),
  };

  for (auto targetSurface : surfaceSequence) {
    step(0.5);
    status();
    BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
    BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
    if (targetSurface.boundary() == 0) {
      BOOST_CHECK_EQUAL(targetSurface.volume(),
                        state.navigation.currentVolume->geometryId().volume());
    }
    target();

    step(1.0);
    status();
    BOOST_CHECK_NE(state.navigation.currentSurface, nullptr);
    if (!state.navigation.navigationBreak) {
      BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
    }
    BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                      targetSurface);

    if (targetSurface.boundary() == 0) {
      BOOST_CHECK_EQUAL(targetSurface.volume(),
                        state.navigation.currentVolume->geometryId().volume());
    }
    target();
  }

  BOOST_CHECK_EQUAL(state.navigation.navigationBreak, true);
  BOOST_CHECK_EQUAL(state.navigation.currentVolume, nullptr);
  BOOST_CHECK_EQUAL(state.navigation.currentSurface->geometryId(),
                    surfaceSequence.back());
}

}  // namespace Test
}  // namespace Acts
