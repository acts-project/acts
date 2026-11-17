// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <utility>

namespace ActsTests {
namespace {
using namespace Acts;

class CountingField final : public MagneticFieldProvider {
 public:
  mutable unsigned int queries = 0;
  bool fail = false;

  Cache makeCache(const MagneticFieldContext& /*context*/) const override {
    return Cache(std::in_place_type<int>, 0);
  }

  Result<Vector3> getField(const Vector3& /*position*/,
                           Cache& /*cache*/) const override {
    ++queries;
    if (fail) {
      return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
    }
    return Result<Vector3>::success(Vector3(0., 0., 2.));
  }
};

struct CustomStepper {
  struct State {
    unsigned int queries = 0;
    unsigned int statusUpdates = 0;
    double curvature = 1.;
    Vector3 usedCurvature = Vector3::Zero();
    bool fail = false;
  };

  Result<Vector3> curvature(State& state) const {
    ++state.queries;
    if (state.fail) {
      return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
    }
    return Result<Vector3>::success(Vector3(state.curvature, 0., 0.));
  }

  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& /*surface*/, std::uint8_t /*index*/,
      Direction /*direction*/, const BoundaryTolerance& /*boundaryTolerance*/,
      double /*surfaceTolerance*/, ConstrainedStep::Type /*stype*/,
      const Logger& /*logger*/, const Vector3& curvature) const {
    ++state.statusUpdates;
    state.usedCurvature = curvature;
    state.curvature *= 2.;
    return IntersectionStatus::onSurface;
  }
};

struct CompositeStepper : CustomStepper {
  Result<IntersectionStatus> updateSurfaceStatusWithCurvature(
      State& state, const Surface& /*surface*/, std::uint8_t /*index*/,
      Direction /*direction*/, const BoundaryTolerance& /*boundaryTolerance*/,
      double /*surfaceTolerance*/, ConstrainedStep::Type /*stype*/,
      const Logger& /*logger*/) const {
    ++state.statusUpdates;
    // The override owns its curvature queries and may change the components.
    const auto value = curvature(state);
    if (!value.ok()) {
      return value.error();
    }
    state.usedCurvature = *value;
    state.curvature *= 3.;
    return Result<IntersectionStatus>::success(IntersectionStatus::onSurface);
  }
};

struct LegacyStepper : CustomStepper {
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& /*surface*/, std::uint8_t /*index*/,
      Direction /*direction*/, const BoundaryTolerance& /*boundaryTolerance*/,
      double /*surfaceTolerance*/, ConstrainedStep::Type /*stype*/,
      const Logger& /*logger*/) const {
    ++state.statusUpdates;
    state.curvature *= 4.;
    return IntersectionStatus::onSurface;
  }
};

static_assert(detail::kSurfaceStatusPreservesKinematics<EigenStepper<>>);
static_assert(!detail::kSurfaceStatusPreservesKinematics<CustomStepper>);
static_assert(!detail::kSurfaceStatusPreservesKinematics<CompositeStepper>);
static_assert(!detail::kSurfaceStatusPreservesKinematics<LegacyStepper>);

struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;
  std::shared_ptr<CountingField> field = std::make_shared<CountingField>();
  EigenStepper<> stepper{field};
  EigenStepper<>::State state =
      stepper.makeState(EigenStepper<>::Options(gctx, mctx));
  std::shared_ptr<PlaneSurface> surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());

  Fixture() {
    const auto parameters = BoundTrackParameters::createCurvilinear(
        Vector4(0., 0., -2., 0.), Vector3(1., 0., 1.).normalized(), 1.,
        std::nullopt, ParticleHypothesis::pion());
    stepper.initialize(state, parameters);
  }

  template <typename Stepper>
  Result<IntersectionStatus> update(
      const Stepper& targetStepper, typename Stepper::State& targetState,
      const Vector3* knownCurvature = nullptr,
      std::optional<Vector3>* evaluatedCurvature = nullptr) const {
    return detail::updateSurfaceStatusWithCurvature(
        targetStepper, targetState, *surface, 0, Direction::Forward(),
        BoundaryTolerance::None(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator, getDummyLogger(), knownCurvature,
        evaluatedCurvature);
  }
};
}  // namespace

BOOST_FIXTURE_TEST_SUITE(SteppingHelper, Fixture)

BOOST_AUTO_TEST_CASE(ReuseEigenStepperCurvature) {
  const auto known = detail::getCurvature(stepper, state);
  BOOST_REQUIRE(known.ok());
  BOOST_CHECK_EQUAL(field->queries, 1u);
  const FreeVector before = state.pars;
  std::optional<Vector3> captured;

  const auto status = update(stepper, state, &*known, &captured);
  BOOST_REQUIRE(status.ok());
  BOOST_CHECK(*status == IntersectionStatus::reachable);
  BOOST_CHECK_EQUAL(field->queries, 1u);
  BOOST_REQUIRE(captured.has_value());
  BOOST_CHECK_EQUAL(*captured, *known);
  BOOST_CHECK_EQUAL(state.pars, before);

  // Reusing the output as the next input must survive clearing the output.
  BOOST_REQUIRE(update(stepper, state, &*captured, &captured).ok());
  BOOST_REQUIRE(captured.has_value());
  BOOST_CHECK_EQUAL(*captured, *known);
  BOOST_CHECK_EQUAL(field->queries, 1u);
}

BOOST_AUTO_TEST_CASE(CaptureFreshEigenStepperCurvatureAtArrival) {
  state.pars.segment<3>(eFreePos0).setZero();
  std::optional<Vector3> captured = Vector3::Constant(99.);
  const FreeVector before = state.pars;

  const auto status = update(stepper, state, nullptr, &captured);
  BOOST_REQUIRE(status.ok());
  BOOST_CHECK(*status == IntersectionStatus::onSurface);
  BOOST_CHECK_EQUAL(field->queries, 1u);
  BOOST_REQUIRE(captured.has_value());
  const Vector3 expected = stepper.qOverP(state) *
                           stepper.direction(state).cross(Vector3(0., 0., 2.));
  BOOST_CHECK_EQUAL(*captured, expected);
  BOOST_CHECK_EQUAL(state.pars, before);

  // A changed state needs a fresh query; the old output is not an input.
  state.pars[eFreeQOverP] *= 2.;
  BOOST_REQUIRE(update(stepper, state, nullptr, &captured).ok());
  BOOST_CHECK_EQUAL(field->queries, 2u);
  BOOST_REQUIRE(captured.has_value());
  BOOST_CHECK_EQUAL(*captured, 2. * expected);
}

BOOST_AUTO_TEST_CASE(FieldFailureClearsCapturedCurvature) {
  field->fail = true;
  std::optional<Vector3> captured = Vector3::Ones();
  const auto status = update(stepper, state, nullptr, &captured);
  BOOST_REQUIRE(!status.ok());
  BOOST_CHECK(status.error() == MagneticFieldError::OutOfBounds);
  BOOST_CHECK_EQUAL(field->queries, 1u);
  BOOST_CHECK(!captured.has_value());
}

BOOST_AUTO_TEST_CASE(CustomStepperKeepsFreshQueries) {
  CustomStepper custom;
  CustomStepper::State customState;
  const Vector3 stale(99., 0., 0.);
  std::optional<Vector3> captured = stale;
  const auto status = update(custom, customState, &stale, &captured);
  BOOST_REQUIRE(status.ok());
  BOOST_CHECK(*status == IntersectionStatus::onSurface);
  BOOST_CHECK_EQUAL(customState.queries, 1u);
  BOOST_CHECK_EQUAL(customState.statusUpdates, 1u);
  BOOST_CHECK_EQUAL(customState.usedCurvature, Vector3(1., 0., 0.));
  BOOST_CHECK(!captured.has_value());
  // Its status update changed kinematics: the pre-update curvature is stale.
  const auto current = detail::getCurvature(custom, customState);
  BOOST_REQUIRE(current.ok());
  BOOST_CHECK_EQUAL(*current, Vector3(2., 0., 0.));

  customState.fail = true;
  captured = stale;
  const auto failure = update(custom, customState, &stale, &captured);
  BOOST_REQUIRE(!failure.ok());
  BOOST_CHECK(failure.error() == MagneticFieldError::OutOfBounds);
  BOOST_CHECK_EQUAL(customState.statusUpdates, 1u);
  BOOST_CHECK(!captured.has_value());
}

BOOST_AUTO_TEST_CASE(CompositeStepperOwnsItsQueries) {
  CompositeStepper composite;
  CompositeStepper::State compositeState;
  const Vector3 stale(99., 0., 0.);
  std::optional<Vector3> captured = stale;
  const auto status = update(composite, compositeState, &stale, &captured);
  BOOST_REQUIRE(status.ok());
  BOOST_CHECK(*status == IntersectionStatus::onSurface);
  BOOST_CHECK_EQUAL(compositeState.queries, 1u);
  BOOST_CHECK_EQUAL(compositeState.statusUpdates, 1u);
  BOOST_CHECK_EQUAL(compositeState.usedCurvature, Vector3(1., 0., 0.));
  BOOST_CHECK_EQUAL(compositeState.curvature, 3.);
  BOOST_CHECK(!captured.has_value());

  compositeState.fail = true;
  captured = stale;
  const auto failure = update(composite, compositeState, &stale, &captured);
  BOOST_REQUIRE(!failure.ok());
  BOOST_CHECK(failure.error() == MagneticFieldError::OutOfBounds);
  BOOST_CHECK(!captured.has_value());
}

BOOST_AUTO_TEST_CASE(LegacyStepperDoesNotGainQueries) {
  LegacyStepper legacy;
  LegacyStepper::State legacyState;
  // Even a failing accessor must not be called on the legacy path.
  legacyState.fail = true;
  const Vector3 stale(99., 0., 0.);
  std::optional<Vector3> captured = stale;
  const auto status = update(legacy, legacyState, &stale, &captured);
  BOOST_REQUIRE(status.ok());
  BOOST_CHECK(*status == IntersectionStatus::onSurface);
  BOOST_CHECK_EQUAL(legacyState.queries, 0u);
  BOOST_CHECK_EQUAL(legacyState.statusUpdates, 1u);
  BOOST_CHECK_EQUAL(legacyState.curvature, 4.);
  BOOST_CHECK(!captured.has_value());
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
