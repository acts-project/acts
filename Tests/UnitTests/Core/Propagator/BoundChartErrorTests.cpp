// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/BoundParametrizationError.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/HelixStepper.hpp"
#include "Acts/Propagator/MultiStepperLoop.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/RiddersStepper.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

using namespace Acts;

namespace ActsTests {
namespace {
const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;

template <typename T>
void checkSingular(const Result<T>& result) {
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == BoundParametrizationError::SingularDirection);
}

template <typename Stepper>
void checkStepper(const Stepper& stepper) {
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  const BoundVector bound(0., 0., 0.3, 1.2, 0.1, 0.);
  const BoundMatrix covariance = BoundMatrix::Identity();
  for (double sign : {-1., 1.}) {
    typename Stepper::Options options(gctx, mctx);
    auto state = stepper.makeState(options);
    stepper.initialize(state, BoundTrackParameters(surface, bound, covariance,
                                                   ParticleHypothesis::pion()));
    const Vector3 direction(0., 0., sign);
    if constexpr (requires { state.primaryStepperState; }) {
      // Set the propagated primary direction exactly, including the south pole.
      state.primaryStepperState.pars.template segment<3>(eFreeDir0) = direction;
    } else {
      stepper.update(state, Vector3::Zero(), direction, 0.1, 0.);
    }
    const auto checkUnchanged = [&] {
      BOOST_REQUIRE(state.cov.has_value());
      BOOST_CHECK_EQUAL(*state.cov, covariance);
      BOOST_CHECK_EQUAL(stepper.position(state), Vector3::Zero());
      BOOST_CHECK_EQUAL(stepper.direction(state), direction);
    };
    checkSingular(stepper.transportToBound(state, *surface));
    checkUnchanged();
    checkSingular(stepper.transportToCurvilinear(state));
    checkUnchanged();
    state.cov.reset();
    BOOST_CHECK(stepper.transportToBound(state, *surface).ok());
    BOOST_CHECK(stepper.transportToCurvilinear(state).ok());
  }
}
}  // namespace

BOOST_AUTO_TEST_CASE(SurfaceChartErrors) {
  const auto plane = Surface::makeShared<PlaneSurface>(Transform3::Identity());
  const auto disc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 1., 5.);
  for (double sign : {-1., 1.}) {
    const Vector3 direction(0., 0., sign);
    const Vector3 position(2., 0., 0.);
    checkSingular(plane->freeToBoundJacobian(gctx, position, direction));
    checkSingular(disc->freeToBoundJacobian(gctx, position, direction));
    checkSingular(CurvilinearSurface(direction).freeToBoundJacobian());
    checkSingular(plane->alignmentToBoundDerivative(gctx, position, direction,
                                                    FreeVector::Zero()));
  }
}

BOOST_AUTO_TEST_CASE(CovarianceFailurePreservesState) {
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  for (double sign : {-1., 1.}) {
    const Vector3 direction(0., 0., sign);
    FreeVector parameters = FreeVector::Zero();
    parameters.segment<3>(eFreeDir0) = direction;
    parameters[eFreeQOverP] = 0.1;
    const FreeVector originalParameters = parameters;
    BoundMatrix covariance = 2. * BoundMatrix::Identity();
    BoundMatrix jacobian = 3. * BoundMatrix::Identity();
    FreeMatrix transport = 4. * FreeMatrix::Identity();
    FreeVector derivatives = FreeVector::Ones();
    BoundToFreeMatrix anchor = BoundToFreeMatrix::Constant(5.);
    const auto checkUnchanged = [&] {
      BOOST_CHECK_EQUAL(covariance, 2. * BoundMatrix::Identity());
      BOOST_CHECK_EQUAL(jacobian, 3. * BoundMatrix::Identity());
      BOOST_CHECK_EQUAL(transport, 4. * FreeMatrix::Identity());
      BOOST_CHECK_EQUAL(derivatives, FreeVector::Ones());
      BOOST_CHECK_EQUAL(anchor, BoundToFreeMatrix::Constant(5.));
      BOOST_CHECK_EQUAL(parameters, originalParameters);
    };
    checkSingular(detail::transportCovarianceToBound(
        gctx, *surface, covariance, jacobian, transport, derivatives, anchor,
        std::nullopt, parameters, FreeToBoundCorrection{}));
    checkUnchanged();
    checkSingular(detail::transportCovarianceToCurvilinear(
        covariance, jacobian, transport, derivatives, anchor, std::nullopt,
        direction));
    checkUnchanged();
    checkSingular(detail::sympy::transportCovarianceToBound(
        gctx, *surface, covariance, jacobian, derivatives, anchor, std::nullopt,
        parameters, FreeToBoundCorrection{}));
    checkUnchanged();
    checkSingular(detail::sympy::transportCovarianceToCurvilinear(
        covariance, jacobian, derivatives, anchor, std::nullopt, direction));
    checkUnchanged();
  }
}

BOOST_AUTO_TEST_CASE(StepperChartErrors) {
  const auto field = std::make_shared<ConstantBField>(Vector3::Zero());
  checkStepper(StraightLineStepper{});
  checkStepper(EigenStepper<>(field));
  checkStepper(HelixStepper(field));
  checkStepper(SympyStepper(field));
  checkStepper(AtlasStepper(field));
  checkStepper(Experimental::RiddersStepper<EigenStepper<>>(field));
}

BOOST_AUTO_TEST_CASE(SympyFailurePreservesMaterialAndScaling) {
  SympyStepper stepper(std::make_shared<ConstantBField>(Vector3::Zero()));
  SympyStepper::Options options(gctx, mctx);
  auto state = stepper.makeState(options);
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  stepper.initialize(state,
                     BoundTrackParameters(
                         surface, BoundVector(0., 0., 0.3, 1.2, 0.1, 0.),
                         BoundMatrix::Identity(), ParticleHypothesis::pion()));
  stepper.update(state, Vector3::Zero(), Vector3::UnitZ(), 0.1, 0.);
  state.materialEffectsAccumulator.initialize(0.01, ParticleHypothesis::pion(),
                                              10.);
  state.materialEffectsAccumulator.accumulate(makeSilicon(), 1., 0.1, 0.1);
  const auto before =
      state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
          Vector3::UnitZ());
  BOOST_REQUIRE(before.has_value());
  const auto anchor = state.jacToGlobal;
  const auto derivatives = state.derivative;
  for (bool curvilinear : {false, true}) {
    checkSingular(curvilinear ? stepper.transportToCurvilinear(state)
                              : stepper.transportToBound(state, *surface));
    BOOST_CHECK_EQUAL(state.jacToGlobal, anchor);
    BOOST_CHECK_EQUAL(state.derivative, derivatives);
    const auto after =
        state.materialEffectsAccumulator.computeAdditionalFreeCovariance(
            Vector3::UnitZ());
    BOOST_REQUIRE(after.has_value());
    BOOST_CHECK_EQUAL(*after, *before);
  }
}

BOOST_AUTO_TEST_CASE(MultiComponentChartErrors) {
  const auto field = std::make_shared<ConstantBField>(Vector3::Zero());
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  const auto check = [&](const auto& stepper) {
    using Stepper = std::decay_t<decltype(stepper)>;
    typename Stepper::Options options(gctx, mctx);
    MultiComponentBoundTrackParameters parameters(surface, true,
                                                  ParticleHypothesis::pion());
    const BoundVector bound(0., 0., 0.3, 1.2, 0.1, 0.);
    const BoundMatrix covariance = BoundMatrix::Identity();
    parameters.pushComponent(0.9, bound, covariance);
    parameters.pushComponent(0.1, bound, covariance);
    for (bool curvilinear : {false, true}) {
      auto state = stepper.makeState(options);
      stepper.initialize(state, parameters);
      // The lower-weight component is deliberately not the reducer's selection.
      state.components[1].state.pars.template segment<3>(eFreeDir0) =
          Vector3::UnitZ();
      checkSingular(curvilinear ? stepper.transportToCurvilinear(state)
                                : stepper.transportToBound(state, *surface));
    }
  };
  check(MultiStepperLoop<EigenStepper<>>(field));
  check(MultiStepperLoop<SympyStepper>(field));
}

BOOST_AUTO_TEST_CASE(PropagatorReportsChartError) {
  Propagator propagator(StraightLineStepper{}, VoidNavigator{});
  decltype(propagator)::Options<> options(gctx, mctx);
  options.pathLimit = 1.;
  for (bool withCovariance : {false, true}) {
    const auto start = BoundTrackParameters::createCurvilinear(
        Vector4::Zero(), 0., 0., 0.1,
        withCovariance ? std::optional<BoundMatrix>(BoundMatrix::Identity())
                       : std::nullopt,
        ParticleHypothesis::pion());
    const auto result = propagator.propagate(start, options);
    if (withCovariance) {
      checkSingular(result);
    } else {
      BOOST_CHECK(result.ok());
    }
  }
}

}  // namespace ActsTests
