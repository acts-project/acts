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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"

#include <cmath>
#include <memory>
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
  bool onPortal = false;
};

struct RecordTransitions {
  struct result_type {
    std::vector<Transition> transitions;
    std::size_t steps = 0;
    std::size_t backtracks = 0;
  };

  template <typename State, typename Stepper, typename Nav>
  Result<void> act(const State& state, const Stepper& stepper, const Nav&,
                   result_type& result, const Logger&) const {
    result.steps = state.steps;
    result.backtracks = state.navigation.statistics.nPortalBacktracks;
    const auto* volume = state.navigation.currentVolume;
    const std::string name = volume == nullptr ? "NONE" : volume->volumeName();
    if (!result.transitions.empty() &&
        result.transitions.back().volume == name) {
      return Result<void>::success();
    }
    const auto* surface = state.navigation.currentSurface;
    result.transitions.push_back(
        {name, state.options.direction * stepper.pathLength(state.stepping),
         stepper.position(state.stepping),
         surface != nullptr && surface->geometryId().boundary() != 0});
    return Result<void>::success();
  }
};

// Two shells stacked along z
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

// A cylinder with a shell around it
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

// Two shells stacked along r, followed along z by a shell that spans both. The
// disc between them is a single portal.
std::shared_ptr<const TrackingGeometry> mergedDisc() {
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {50_mm, 50_mm};
  config.envelope[AxisDirection::AxisZ] = {100_mm, 100_mm};
  Blueprint root(config);
  root.addCylinderContainer("Stack", AxisDirection::AxisZ, [](auto& stack) {
    stack.addCylinderContainer("Radial", AxisDirection::AxisR, [](auto& rad) {
      rad.addStaticVolume(
          Transform3{Translation3{0., 0., -50_mm}},
          std::make_shared<CylinderVolumeBounds>(50_mm, 100_mm, 50_mm), "A");
      rad.addStaticVolume(
          Transform3{Translation3{0., 0., -50_mm}},
          std::make_shared<CylinderVolumeBounds>(100_mm, 150_mm, 50_mm), "B");
    });
    stack.addStaticVolume(
        Transform3{Translation3{0., 0., 50_mm}},
        std::make_shared<CylinderVolumeBounds>(50_mm, 150_mm, 50_mm), "C");
  });
  return root.construct({}, gctx);
}

RecordTransitions::result_type propagate(
    const std::shared_ptr<const TrackingGeometry>& geometry,
    const Vector3& position, const Vector3& direction, double qOverP,
    double pathLimit, bool backtracking, double bField = fieldZ,
    Direction propagationDirection = Direction::Forward()) {
  Propagator propagator{
      EigenStepper<>(std::make_shared<ConstantBField>(Vector3{0., 0., bField})),
      Navigator({.trackingGeometry = geometry})};
  decltype(propagator)::Options<ActorList<RecordTransitions, EndOfWorldReached>>
      options(gctx, mctx);
  options.direction = propagationDirection;
  options.pathLimit = pathLimit;
  options.maxSteps = 1000;
  options.stepping.maxStepSize = 100_mm;
  options.navigation.portalBacktracking = backtracking;
  // Reversing both the momentum and the charge gives the same path when
  // propagating backward
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

// Helix in a constant field along z, omega = (q/p) B
Vector3 helixPosition(const Vector3& position, const Vector3& direction,
                      double omega, double path) {
  const double sine = std::sin(omega * path) / omega;
  const double cosine = (1. - std::cos(omega * path)) / omega;
  return position + Vector3{sine * direction.x() + cosine * direction.y(),
                            sine * direction.y() - cosine * direction.x(),
                            path * direction.z()};
}

// Path at which the helix crosses the radius, between two paths that
// bracket it
double radiusCrossing(const Vector3& position, const Vector3& direction,
                      double omega, double radius, double lower, double upper) {
  auto distance = [&](double path) {
    return helixPosition(position, direction, omega, path).head<2>().norm() -
           radius;
  };
  BOOST_REQUIRE(distance(lower) * distance(upper) < 0.);
  for (int i = 0; i < 100; ++i) {
    const double middle = 0.5 * (lower + upper);
    if (distance(lower) * distance(middle) <= 0.) {
      upper = middle;
    } else {
      lower = middle;
    }
  }
  return 0.5 * (lower + upper);
}

void checkTransition(const Transition& actual, const std::string& volume,
                     double path, const Vector3& position) {
  BOOST_CHECK_EQUAL(actual.volume, volume);
  BOOST_CHECK_SMALL(actual.path - path, 1_um);
  BOOST_CHECK_SMALL((actual.position - position).norm(), 1_um);
  BOOST_CHECK(actual.onPortal);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(NavigatorPortalBacktracking)

// The straight line from the start misses the disc, and reaches the outer
// cylinder beyond it instead
BOOST_AUTO_TEST_CASE(DiscMissedByStraightLine) {
  const auto geometry = axialStack();
  const Vector3 start{100_mm, 0., -50_mm};
  const Vector3 direction{0., -0.8, 0.6};
  // A circle of 100 mm around the axis stays within the shell
  const double omega = 0.8 / 100_mm;
  const double discPath = 50_mm / direction.z();

  for (const auto propagationDirection :
       {Direction::Forward(), Direction::Backward()}) {
    const auto result = propagate(geometry, start, direction, omega / fieldZ,
                                  150_mm, true, fieldZ, propagationDirection);
    BOOST_REQUIRE_EQUAL(result.transitions.size(), 2u);
    BOOST_CHECK_EQUAL(result.transitions.front().volume, "first");
    checkTransition(result.transitions[1], "second", discPath,
                    helixPosition(start, direction, omega, discPath));
    BOOST_CHECK_EQUAL(result.backtracks, 1u);
  }
}

// The straight line from the start misses the inner cylinder, which the
// trajectory bends into
BOOST_AUTO_TEST_CASE(InnerCylinderMissedByStraightLine) {
  const auto geometry = radialStack();
  const Vector3 start{125_mm, 0., 0.};
  const Vector3 direction{0., -1., 0.};
  // A circle of 75 mm around (50,0)
  const double omega = 1. / 75_mm;
  const double crossingPath =
      std::acos((60. * 60. - 50. * 50. - 75. * 75.) / (2. * 50. * 75.)) / omega;

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

// The trajectory dips into the inner cylinder and comes back out within a
// single step
BOOST_AUTO_TEST_CASE(DipIntoInnerCylinder) {
  const auto geometry = radialStack();
  // A circle of 200 mm around (141,0) encloses the axis and gets as close as
  // 59 mm to it, 60 mm of path after the start. Its tangents stay further
  // away, so the straight line from the start misses the inner cylinder.
  const double radius = 200_mm;
  const double omega = 1. / radius;
  const double angle = 60_mm / radius;
  const Vector3 start{141_mm - radius * std::cos(angle),
                      -radius * std::sin(angle), 0.};
  const Vector3 direction{-std::sin(angle), std::cos(angle), 0.};
  BOOST_CHECK_CLOSE(helixPosition(start, direction, omega, 60_mm).x(), -59_mm,
                    1e-6);
  const double entryPath =
      radiusCrossing(start, direction, omega, 60_mm, 0_mm, 60_mm);
  const double exitPath =
      radiusCrossing(start, direction, omega, 60_mm, 60_mm, 120_mm);

  const auto result =
      propagate(geometry, start, direction, omega / fieldZ, 150_mm, true);
  BOOST_REQUIRE_EQUAL(result.transitions.size(), 3u);
  BOOST_CHECK_EQUAL(result.transitions[0].volume, "outer");
  checkTransition(result.transitions[1], "inner", entryPath,
                  helixPosition(start, direction, omega, entryPath));
  checkTransition(result.transitions[2], "outer", exitPath,
                  helixPosition(start, direction, omega, exitPath));
}

// The trajectory leaves through the outer cylinder right before a disc that
// extends over the neighbouring volume as well
BOOST_AUTO_TEST_CASE(PortalBeyondTheVolumeFace) {
  const auto geometry = mergedDisc();
  const Vector3 start{90_mm, 0., -50_mm};
  const Vector3 direction{0., -0.6, 0.8};
  // A circle of 200 mm around (290,0), seen along the axis
  const double omega = -0.6 / 200_mm;
  const double cylinderPath =
      radiusCrossing(start, direction, omega, 100_mm, 0_mm, 62.5_mm);
  const double discPath = 50_mm / direction.z();
  BOOST_CHECK_LT(cylinderPath, discPath);

  const auto result =
      propagate(geometry, start, direction, omega / fieldZ, 100_mm, true);
  BOOST_REQUIRE_EQUAL(result.transitions.size(), 3u);
  BOOST_CHECK_EQUAL(result.transitions[0].volume, "A");
  checkTransition(result.transitions[1], "B", cylinderPath,
                  helixPosition(start, direction, omega, cylinderPath));
  checkTransition(result.transitions[2], "C", discPath,
                  helixPosition(start, direction, omega, discPath));
}

// Straight lines do not leave a volume unnoticed
BOOST_AUTO_TEST_CASE(StraightLineIsUnchanged) {
  for (const auto& geometry : {axialStack(), radialStack(), mergedDisc()}) {
    for (const Vector3& direction :
         {Vector3{0., -0.8, 0.6}, Vector3{-0.6, 0., 0.8},
          Vector3{0.6, 0.64, -0.48}}) {
      const Vector3 start{100_mm, 0., -50_mm};
      const auto disabled =
          propagate(geometry, start, direction, 1. / 1_GeV, 500_mm, false, 0.);
      const auto enabled =
          propagate(geometry, start, direction, 1. / 1_GeV, 500_mm, true, 0.);
      BOOST_CHECK_EQUAL(enabled.backtracks, 0u);
      BOOST_CHECK_EQUAL(enabled.steps, disabled.steps);
      BOOST_REQUIRE_EQUAL(enabled.transitions.size(),
                          disabled.transitions.size());
      for (std::size_t i = 0; i < enabled.transitions.size(); ++i) {
        BOOST_CHECK_EQUAL(enabled.transitions[i].volume,
                          disabled.transitions[i].volume);
        BOOST_CHECK_EQUAL(enabled.transitions[i].path,
                          disabled.transitions[i].path);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
