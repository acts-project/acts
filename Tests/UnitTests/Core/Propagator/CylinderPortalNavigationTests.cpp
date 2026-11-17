// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ActsTests {
namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;
constexpr double fieldZ = 2_T;

struct Transition {
  std::string volume;
  double path = 0.;
  Vector3 position = Vector3::Zero();
  bool onFiniteSurface = false;
  bool hasPendingTarget = false;
};

// Observe the ordinary post-step actor interface: a supporting-plane waypoint
// must not masquerade as an actual portal arrival.
struct RecordTransitions {
  struct result_type {
    std::vector<Transition> transitions;
    std::vector<double> surfacePaths;
    std::vector<const Surface*> reachedSurfaces;
    std::size_t steps = 0;
  };

  template <typename State, typename Stepper, typename Nav>
  Result<void> act(const State& state, const Stepper& stepper, const Nav&,
                   result_type& result, const Logger&) const {
    if (state.stage == PropagatorStage::preStep) {
      return Result<void>::success();
    }
    result.steps = state.steps;
    if (state.stage == PropagatorStage::postStep &&
        state.navigation.currentSurface != nullptr) {
      result.surfacePaths.push_back(state.options.direction *
                                    stepper.pathLength(state.stepping));
      result.reachedSurfaces.push_back(state.navigation.currentSurface);
    }
    const auto* volume = state.navigation.currentVolume;
    const std::string name = volume == nullptr ? "NONE" : volume->volumeName();
    if (!result.transitions.empty() &&
        result.transitions.back().volume == name) {
      return Result<void>::success();
    }
    const Vector3 position = stepper.position(state.stepping);
    const auto* surface = state.navigation.currentSurface;
    result.transitions.push_back(
        {name, state.options.direction * stepper.pathLength(state.stepping),
         position,
         surface != nullptr &&
             surface->isOnSurface(gctx, position,
                                  stepper.direction(state.stepping),
                                  BoundaryTolerance::None(), 1_um),
         state.navigation.pendingTarget.has_value()});
    return Result<void>::success();
  }
};

std::shared_ptr<const TrackingGeometry> axialStack() {
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {90_mm, 90_mm};
  config.envelope[AxisDirection::AxisZ] = {100_mm, 100_mm};
  Blueprint root(config);
  root.addCylinderContainer("Stack", AxisDirection::AxisZ, [](auto& stack) {
    stack.addStaticVolume(
        Transform3{Translation3{0., 0., -50_mm}},
        std::make_shared<CylinderVolumeBounds>(90_mm, 110_mm, 50_mm), "first");
    stack.addStaticVolume(
        Transform3{Translation3{0., 0., 50_mm}},
        std::make_shared<CylinderVolumeBounds>(90_mm, 110_mm, 50_mm), "second");
  });
  return root.construct({}, gctx);
}

std::shared_ptr<const TrackingGeometry> radialStack() {
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {0_mm, 50_mm};
  config.envelope[AxisDirection::AxisZ] = {100_mm, 100_mm};
  Blueprint root(config);
  root.addCylinderContainer("Stack", AxisDirection::AxisR, [](auto& stack) {
    stack.addStaticVolume(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0_mm, 60_mm, 500_mm), "inner");
    stack.addStaticVolume(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(60_mm, 150_mm, 500_mm), "outer");
  });
  return root.construct({}, gctx);
}

std::shared_ptr<const TrackingGeometry> finiteChild() {
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {0_mm, 100_mm};
  config.envelope[AxisDirection::AxisZ] = {150_mm, 150_mm};
  Blueprint root(config);
  root.addStaticVolume(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 50_mm), "child");
  return root.construct({}, gctx);
}

RecordTransitions::result_type propagate(
    const std::shared_ptr<const TrackingGeometry>& geometry,
    const Vector3& position, const Vector3& direction, double qOverP,
    double pathLimit, std::optional<bool> enableGuard, double bField = fieldZ,
    Direction propagationDirection = Direction::Forward(),
    const Surface* additionalSurface = nullptr) {
  Propagator propagator{
      EigenStepper<>(std::make_shared<ConstantBField>(Vector3{0., 0., bField})),
      Navigator({.trackingGeometry = geometry})};
  decltype(propagator)::Options<ActorList<RecordTransitions, EndOfWorldReached>>
      options(gctx, mctx);
  options.direction = propagationDirection;
  options.pathLimit = pathLimit;
  options.loopProtection = false;
  options.maxSteps = 1000;
  options.stepping.maxStepSize = 100_mm;
  if (additionalSurface != nullptr) {
    options.navigation.registerAdditionalSurface(*additionalSurface);
  }
  if (enableGuard.has_value()) {
    options.navigation.experimentalCylindricalPortalGuard = *enableGuard;
  }
  // Reversing both initial momentum and charge gives the same geometric path
  // under backward propagation, exercising its signed path/step conventions.
  const auto parameters = BoundTrackParameters::createCurvilinear(
      Vector4{position.x(), position.y(), position.z(), 0.},
      propagationDirection * direction, propagationDirection * qOverP,
      std::nullopt, ParticleHypothesis::muon());
  const auto result = propagator.propagate(parameters, options, false);
  if (!result.ok()) {
    BOOST_FAIL(result.error().message());
  }
  return result->get<RecordTransitions::result_type>();
}

const TrackingVolume& startingVolume(const TrackingGeometry& geometry,
                                     const Vector3& position) {
  const auto result = geometry.resolveLowestTrackingVolume(gctx, position);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(*result != nullptr);
  return **result;
}

void checkTransition(const Transition& actual, const std::string& volume,
                     double path, const Vector3& position) {
  BOOST_CHECK_EQUAL(actual.volume, volume);
  BOOST_CHECK_SMALL(actual.path - path, 1_um);
  BOOST_CHECK_SMALL((actual.position - position).norm(), 1_um);
  BOOST_CHECK(actual.onFiniteSurface);
  BOOST_CHECK(!actual.hasPendingTarget);
}

// Exact solution for constant axial B, independent of the guard's conservative
// quadratic bounds. omega = (q/p) B has inverse-length units in ACTS.
Vector3 helixPosition(const Vector3& position, const Vector3& direction,
                      double omega, double path) {
  const double sine = std::sin(omega * path) / omega;
  const double cosine = (1. - std::cos(omega * path)) / omega;
  return position + Vector3{sine * direction.x() + cosine * direction.y(),
                            sine * direction.y() - cosine * direction.x(),
                            path * direction.z()};
}

void checkIdentical(const RecordTransitions::result_type& lhs,
                    const RecordTransitions::result_type& rhs) {
  BOOST_CHECK_EQUAL(lhs.steps, rhs.steps);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      lhs.surfacePaths.begin(), lhs.surfacePaths.end(),
      rhs.surfacePaths.begin(), rhs.surfacePaths.end());
  BOOST_REQUIRE_EQUAL(lhs.transitions.size(), rhs.transitions.size());
  for (std::size_t i = 0; i < lhs.transitions.size(); ++i) {
    BOOST_CHECK_EQUAL(lhs.transitions[i].volume, rhs.transitions[i].volume);
    BOOST_CHECK_EQUAL(lhs.transitions[i].path, rhs.transitions[i].path);
    BOOST_CHECK(lhs.transitions[i].position == rhs.transitions[i].position);
    BOOST_CHECK_EQUAL(lhs.transitions[i].onFiniteSurface,
                      rhs.transitions[i].onFiniteSurface);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(CylinderPortalNavigation)

BOOST_AUTO_TEST_CASE(FiniteCapRejectedByInitialRay) {
  const auto geometry = axialStack();
  const Vector3 start{100_mm, 0., -50_mm};
  const Vector3 direction{0., -0.8, 0.6};
  // A 100 mm transverse circle centered on the cylinder axis: radius stays
  // inside [90,110] mm, although the initial ray misses the cap's finite
  // bounds.
  const double omega = 0.8 / 100_mm;
  const double capPath = 50_mm / direction.z();
  const auto& volume = startingVolume(*geometry, start);
  BOOST_CHECK_EQUAL(volume.volumeName(), "first");
  bool checkedCap = false;
  for (const auto& portal : volume.portals()) {
    if (portal.surface().type() == Surface::Disc &&
        std::abs(portal.surface().center(gctx).z()) < 1_um) {
      BOOST_CHECK(
          !portal.surface()
               .intersect(gctx, start, direction, BoundaryTolerance::None())
               .closestForward()
               .isValid());
      checkedCap = true;
    }
  }
  BOOST_REQUIRE(checkedCap);

  for (const auto propagationDirection :
       {Direction::Forward(), Direction::Backward()}) {
    const auto result = propagate(geometry, start, direction, omega / fieldZ,
                                  150_mm, true, fieldZ, propagationDirection);
    BOOST_REQUIRE_EQUAL(result.transitions.size(), 2u);
    BOOST_CHECK_EQUAL(result.transitions.front().volume, "first");
    checkTransition(result.transitions[1], "second", capPath,
                    helixPosition(start, direction, omega, capPath));
  }
}

BOOST_AUTO_TEST_CASE(CapPromotionDiscardsStagedAdditionalSurfaceTarget) {
  const auto geometry = axialStack();
  const Vector3 start{100_mm, 0., -50_mm};
  const Vector3 direction{0., -0.8, 0.6};
  const double omega = 0.8 / 100_mm;
  const double capPath = 50_mm / direction.z();

  // The tilted plane is 50 mm away along the initial tangent, closer than
  // the staged outer-cylinder root (~57 mm). The helix bends away from the
  // plane initially and reaches it only after ~88 mm, beyond the omitted cap
  // at 83.33 mm. Its pending geometry target therefore belongs to the old
  // volume when the guard promotes the cap to the active target.
  const Vector3 normal = Vector3{-1.04, 1., 0.}.normalized();
  Transform3 transform = Transform3::Identity();
  transform.linear().col(0) = Vector3::UnitZ();
  transform.linear().col(1) = normal.cross(Vector3::UnitZ());
  transform.linear().col(2) = normal;
  transform.translation() = Vector3{100_mm, -40_mm, 0.};
  const auto additional = Surface::makeShared<PlaneSurface>(
      transform, std::make_shared<RectangleBounds>(500_mm, 500_mm));

  Navigator navigator({.trackingGeometry = geometry});
  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  auto state = navigator.makeState(options);
  BOOST_REQUIRE(
      navigator.initialize(state, {.position = start, .direction = direction})
          .ok());
  const auto target = navigator.nextTarget(state, start, direction);
  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), additional.get());
  BOOST_REQUIRE(state.pendingTarget.has_value());
  BOOST_CHECK(state.pendingTarget->isPortalTarget());

  const auto result =
      propagate(geometry, start, direction, omega / fieldZ, 150_mm, true,
                fieldZ, Direction::Forward(), additional.get());
  BOOST_REQUIRE_EQUAL(result.transitions.size(), 2u);
  checkTransition(result.transitions[1], "second", capPath,
                  helixPosition(start, direction, omega, capPath));

  // Solve 1.04*(1-cos(theta)) - sin(theta) + 0.4 = 0 for its first
  // positive root. The plane must still be reached after the portal switch.
  const double planePath =
      (std::atan2(1., 1.04) - std::acos(1.44 / std::hypot(1.04, 1.))) / omega;
  BOOST_REQUIRE_EQUAL(result.reachedSurfaces.size(), 2u);
  BOOST_CHECK_NE(result.reachedSurfaces[0], additional.get());
  BOOST_CHECK_EQUAL(result.reachedSurfaces[1], additional.get());
  BOOST_CHECK_SMALL(result.surfacePaths[0] - capPath, 1_um);
  BOOST_CHECK_SMALL(result.surfacePaths[1] - planePath, 1_um);
}

BOOST_AUTO_TEST_CASE(InnerCylinderRootAppearsAfterEntry) {
  const auto geometry = radialStack();
  const Vector3 start{125_mm, 0., 0.};
  const Vector3 direction{0., -1., 0.};
  // The transverse circle has center (50,0) and radius 75 mm. Its radius
  // about the geometry axis falls from 125 to 25 mm, crossing the 60 mm face.
  // At entry the tangent ray has no inner-cylinder root at all.
  const double omega = 1. / 75_mm;
  const double crossingPath =
      std::acos((60. * 60. - 50. * 50. - 75. * 75.) / (2. * 50. * 75.)) / omega;
  const auto& volume = startingVolume(*geometry, start);
  bool checkedInner = false;
  for (const auto& portal : volume.portals()) {
    const auto* bounds =
        dynamic_cast<const CylinderBounds*>(&portal.surface().bounds());
    if (bounds != nullptr && bounds->get(CylinderBounds::eR) == 60_mm) {
      BOOST_CHECK(
          !portal.surface()
               .intersect(gctx, start, direction, BoundaryTolerance::None())
               .closestForward()
               .isValid());
      checkedInner = true;
    }
  }
  BOOST_REQUIRE(checkedInner);
  for (const auto propagationDirection :
       {Direction::Forward(), Direction::Backward()}) {
    const auto result = propagate(geometry, start, direction, omega / fieldZ,
                                  200_mm, true, fieldZ, propagationDirection);
    BOOST_REQUIRE_EQUAL(result.transitions.size(), 2u);
    BOOST_CHECK_EQUAL(result.transitions.front().volume, "outer");
    checkTransition(result.transitions[1], "inner", crossingPath,
                    helixPosition(start, direction, omega, crossingPath));
  }
}

BOOST_AUTO_TEST_CASE(ChildCylinderOutsideItsFiniteZExtent) {
  const auto geometry = finiteChild();
  const Vector3 start{150_mm, 0., -150_mm};
  const Vector3 direction{-0.6, 0., 0.8};
  const double qOverP = 1. / 0.3_GeV;
  // The trajectory crosses the child's supporting cylinder below z=-50 mm.
  // This must not produce a portal hit. Entry and exit are through its caps.
  const auto result =
      propagate(geometry, start, direction, qOverP, 300_mm, true);
  BOOST_REQUIRE_EQUAL(result.transitions.size(), 3u);
  BOOST_REQUIRE_EQUAL(result.surfacePaths.size(), 2u);
  BOOST_CHECK_SMALL(result.surfacePaths[0] - 125_mm, 1_um);
  BOOST_CHECK_SMALL(result.surfacePaths[1] - 250_mm, 1_um);
  BOOST_CHECK_EQUAL(result.transitions.front().volume,
                    geometry->highestTrackingVolume()->volumeName());
  checkTransition(result.transitions[1], "child", 125_mm,
                  helixPosition(start, direction, qOverP * fieldZ, 125_mm));
  checkTransition(result.transitions[2],
                  geometry->highestTrackingVolume()->volumeName(), 250_mm,
                  helixPosition(start, direction, qOverP * fieldZ, 250_mm));
}

BOOST_AUTO_TEST_CASE(ZeroFieldPreservesStepsAndArrivals) {
  const auto geometry = finiteChild();
  const Vector3 start{150_mm, 0., -150_mm};
  const Vector3 direction{-0.6, 0., 0.8};
  for (const auto propagationDirection :
       {Direction::Forward(), Direction::Backward()}) {
    const auto disabled = propagate(geometry, start, direction, 1. / 1_GeV,
                                    300_mm, false, 0., propagationDirection);
    const auto enabled = propagate(geometry, start, direction, 1. / 1_GeV,
                                   300_mm, true, 0., propagationDirection);
    BOOST_REQUIRE_EQUAL(disabled.transitions.size(), 3u);
    checkIdentical(disabled, enabled);
  }
}

BOOST_AUTO_TEST_CASE(GuardIsDisabledByDefault) {
  const auto geometry = axialStack();
  const Vector3 start{100_mm, 0., -50_mm};
  const Vector3 direction{0., -0.8, 0.6};
  const double qOverP = (0.8 / 100_mm) / fieldZ;
  const auto defaults =
      propagate(geometry, start, direction, qOverP, 150_mm, std::nullopt);
  const auto disabled =
      propagate(geometry, start, direction, qOverP, 150_mm, false);
  checkIdentical(defaults, disabled);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
