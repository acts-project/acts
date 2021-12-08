// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"

using namespace Acts;
using namespace Acts::VectorHelpers;

MagneticFieldContext magCtx;
GeometryContext geoCtx;

using MultiStepper =
    MultiEigenStepperLoop<StepperExtensionList<DefaultExtension>>;
using SingleStepper = EigenStepper<StepperExtensionList<DefaultExtension>>;

struct Options {
  double tolerance = 1e-4;
  double stepSizeCutOff = 0.0;
  std::size_t maxRungeKuttaStepTrials = 10;
  double mass = 1.0;
};

struct Navigation {};

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_charged_no_cov) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, BoundVector::Ones(), std::nullopt});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  // Test charged parameters without covariance matrix
  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  MultiStepper::State state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                            tolerance);

  MultiStepper ms((bField));

  BOOST_CHECK_EQUAL(cmps.size(), ms.numberComponents(state));

  // Test the result & compare with the input/test for reasonable members
  for (const auto &cmp : state.components) {
    BOOST_CHECK_EQUAL(cmp.state.jacToGlobal, BoundToFreeMatrix::Zero());
    BOOST_CHECK_EQUAL(cmp.state.jacTransport, FreeMatrix::Identity());
    BOOST_CHECK_EQUAL(cmp.state.derivative, FreeVector::Zero());
    BOOST_CHECK_EQUAL(cmp.state.cov, BoundSymMatrix::Zero());
  }

  BOOST_CHECK(not state.covTransport);
  for (const auto &cmp : state.components) {
    BOOST_CHECK(not cmp.state.covTransport);
  }

  BOOST_CHECK_EQUAL(state.navDir, ndir);
  for (const auto &cmp : state.components) {
    BOOST_CHECK_EQUAL(cmp.state.navDir, ndir);
  }

  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  for (const auto &cmp : state.components) {
    BOOST_CHECK_EQUAL(cmp.state.pathAccumulated, 0.);
  }

  BOOST_CHECK_EQUAL(ms.charge(state), 1.);
  for (const auto &cmp : state.components) {
    BOOST_CHECK_EQUAL(cmp.state.q, 1.);
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_neutral_no_cov) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, BoundVector::Ones(), std::nullopt});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  // Test neutral parameters without covariance matrix
  MultiComponentBoundTrackParameters<Neutral> multi_pars(surface, cmps);

  MultiStepper::State state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                            tolerance);

  MultiStepper ms((bField));

  BOOST_CHECK_EQUAL(ms.charge(state), 0.);
  for (const auto &cmp : state.components) {
    BOOST_CHECK_EQUAL(cmp.state.q, 0.);
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_charged_with_cov) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, BoundVector::Ones(), BoundSymMatrix::Identity()});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  // Test charged parameters without covariance matrix
  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  MultiStepper::State state(geoCtx, magCtx, bField, multi_pars, ndir, stepSize,
                            tolerance);

  MultiStepper ms((bField));

  // Test the result & compare with the input/test for reasonable members
  for (const auto &cmp : state.components) {
    BOOST_CHECK_NE(cmp.state.cov, BoundSymMatrix::Zero());
  }

  BOOST_CHECK(state.covTransport);
  for (const auto &cmp : state.components) {
    BOOST_CHECK(cmp.state.covTransport);
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_stepper_state_invalid) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  // Empty component vector
  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  BOOST_CHECK_THROW(MultiStepper::State(geoCtx, magCtx, bField, multi_pars,
                                        ndir, stepSize, tolerance),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(multi_eigen_vs_single_eigen) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;
  const BoundVector pars = BoundVector::Ones();
  const BoundSymMatrix cov = BoundSymMatrix::Identity();

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  SingleBoundTrackParameters<SinglyCharged> single_pars(surface, pars, cov);

  MultiStepper::State multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                  stepSize, tolerance);
  SingleStepper::State single_state(geoCtx, bField->makeCache(magCtx),
                                    single_pars, ndir, stepSize, tolerance);

  MultiStepper multi_stepper(bField);
  SingleStepper single_stepper(bField);

  for (int i = 0; i < 10; ++i) {
    // Single stepper
    struct SinglePropState {
      SingleStepper::State &stepping;
      Options options;
    } single_prop_state{single_state, Options{}};

    auto single_result = single_stepper.step(single_prop_state);
    single_stepper.transportCovarianceToCurvilinear(single_state);

    // Multi stepper;
    struct MultiPropState {
      MultiStepper::State &stepping;
      Options options;
      Navigation navigation;
      GeometryContext geoContext;
    } multi_prop_state{multi_state, Options{}, Navigation{}, geoCtx};

    auto multi_result = multi_stepper.step(multi_prop_state);
    multi_stepper.transportCovarianceToCurvilinear(multi_state);

    // Check equality
    BOOST_REQUIRE(multi_result.ok() == true);
    BOOST_REQUIRE(multi_result.ok() == single_result.ok());

    BOOST_CHECK_EQUAL(*single_result, *multi_result);

    for (const auto &cmp : multi_state.components) {
      BOOST_CHECK_EQUAL(cmp.state.pars, single_state.pars);
      BOOST_CHECK_EQUAL(cmp.state.q, single_state.q);
      BOOST_CHECK_EQUAL(cmp.state.cov, single_state.cov);
      BOOST_CHECK_EQUAL(cmp.state.jacTransport, single_state.jacTransport);
      BOOST_CHECK_EQUAL(cmp.state.jacToGlobal, single_state.jacToGlobal);
      BOOST_CHECK_EQUAL(cmp.state.derivative, single_state.derivative);
    }
  }
}

BOOST_AUTO_TEST_CASE(multi_eigen_component_iterable) {
  const NavigationDirection ndir = backward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps;
  for (int i = 0; i < 4; ++i) {
    cmps.push_back({0.25, BoundVector::Random(), BoundSymMatrix::Random()});
  }

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3::Ones().normalized());

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);

  MultiStepper::State multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                  stepSize, tolerance);

  MultiStepper multi_stepper(bField);
  const auto const_iterable = multi_stepper.constComponentIterable(multi_state);

  BOOST_CHECK_EQUAL(cmps.size(), multi_stepper.numberComponents(multi_state));

  {
    auto proxy_it = const_iterable.begin();
    auto vec_it = cmps.begin();
    for (; vec_it != cmps.end(); ++proxy_it, ++vec_it) {
      BOOST_CHECK((*proxy_it).weight() == std::get<double>(*vec_it));
    }
  }

  for (auto cmp : multi_stepper.componentIterable(multi_state)) {
    cmp.weight() *= 2.0;
    cmp.pars() *= 2.0;
    cmp.cov() *= 2.0;
  }

  {
    auto proxy_it = const_iterable.begin();
    auto i = 0ul;
    for (; proxy_it != const_iterable.end(); ++proxy_it, ++i) {
      const auto [weight, single_pars_for_check] = multi_pars[i];
      BOOST_CHECK_CLOSE((*proxy_it).weight(), 2.0 * weight, 1.e-10);
      BOOST_CHECK((*proxy_it).pars().template segment<4>(eFreePos0).isApprox(
          2.0 * single_pars_for_check.fourPosition(geoCtx), 1.e-10));
      BOOST_CHECK((*proxy_it).pars().template segment<3>(eFreeDir0).isApprox(
          2.0 * single_pars_for_check.unitDirection(), 1.e-10));
      BOOST_CHECK_CLOSE((*proxy_it).pars()[eFreeQOverP],
                        2.0 * single_pars_for_check.charge() /
                            single_pars_for_check.absoluteMomentum(),
                        1.e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_surface_status) {
  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto start_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});

  auto right_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3{1.0, 0.0, 0.0}, Vector3{1.0, 0.0, 0.0});

  auto bField = std::make_shared<NullBField>();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(2, {0.5, BoundVector::Zero(), std::nullopt});
  std::get<BoundVector>(cmps[0])[eBoundTheta] = M_PI_2;
  std::get<BoundVector>(cmps[1])[eBoundTheta] = -M_PI_2;
  std::get<BoundVector>(cmps[0])[eBoundQOverP] = 1.0;
  std::get<BoundVector>(cmps[1])[eBoundQOverP] = 1.0;

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(start_surface,
                                                               cmps);

  BOOST_REQUIRE(std::get<1>(multi_pars[0])
                    .unitDirection()
                    .isApprox(Vector3{1.0, 0.0, 0.0}, 1.e-10));
  BOOST_REQUIRE(std::get<1>(multi_pars[1])
                    .unitDirection()
                    .isApprox(Vector3{-1.0, 0.0, 0.0}, 1.e-10));

  MultiStepper::State multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                  stepSize, tolerance);
  MultiStepper multi_stepper(bField);

  // Update surface status and check
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *right_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::reachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::reachable);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::missed);
  }

  // Step forward now
  {
    struct MultiPropState {
      MultiStepper::State &stepping;
      Options options;
      Navigation navigation;
      GeometryContext geoContext;
    } multi_prop_state{multi_state, Options{}, Navigation{}, geoCtx};

    multi_stepper.step(multi_prop_state);
  }

  // Update surface status and check again
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *right_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::onSurface);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::onSurface);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::missed);
  }

  // Start surface should be unreachable
  {
    auto status =
        multi_stepper.updateSurfaceStatus(multi_state, *start_surface, false);

    BOOST_CHECK(status == Intersection3D::Status::unreachable);

    auto cmp_iterable = multi_stepper.constComponentIterable(multi_state);

    BOOST_CHECK((*cmp_iterable.begin()).status() ==
                Intersection3D::Status::unreachable);
    BOOST_CHECK((*(++cmp_iterable.begin())).status() ==
                Intersection3D::Status::unreachable);
  }
}

BOOST_AUTO_TEST_CASE(test_bound_state) {
  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});
  auto bField = std::make_shared<NullBField>();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSymMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  MultiStepper::State multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                  stepSize, tolerance);
  MultiStepper multi_stepper(bField);

  auto res = multi_stepper.boundState(multi_state, *surface, true);

  BOOST_REQUIRE(res.ok());

  const auto [bound_pars, jacobian, pathLength] = *res;

  BOOST_CHECK(jacobian == decltype(jacobian)::Zero());
  BOOST_CHECK(pathLength == 0.0);
  BOOST_CHECK(bound_pars.parameters().isApprox(pars, 1.e-8));
  BOOST_CHECK(bound_pars.covariance()->isApprox(cov, 1.e-8));
}

BOOST_AUTO_TEST_CASE(test_curvilinear_state) {
  const NavigationDirection ndir = forward;
  const double stepSize = 123.;
  const double tolerance = 234.;

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Vector3::Zero(), Vector3{1.0, 0.0, 0.0});
  auto bField = std::make_shared<NullBField>();

  // Use Ones() here, so that the angles are in correct range
  const auto pars = BoundVector::Ones().eval();
  const auto cov = []() {
    auto c = BoundSymMatrix::Random().eval();
    c *= c.transpose();
    return c;
  }();

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSymMatrix>>>
      cmps(4, {0.25, pars, cov});
  SingleBoundTrackParameters<SinglyCharged> check_pars(surface, pars, cov);

  MultiComponentBoundTrackParameters<SinglyCharged> multi_pars(surface, cmps);
  MultiStepper::State multi_state(geoCtx, magCtx, bField, multi_pars, ndir,
                                  stepSize, tolerance);
  MultiStepper multi_stepper(bField);

  const auto [curv_pars, jac, pathLength] =
      multi_stepper.curvilinearState(multi_state);

  BOOST_CHECK(
      curv_pars.fourPosition(multi_state.geoContext)
          .isApprox(check_pars.fourPosition(multi_state.geoContext), 1.e-8));
  BOOST_CHECK(
      curv_pars.unitDirection().isApprox(check_pars.unitDirection(), 1.e-8));
  BOOST_CHECK_CLOSE(curv_pars.absoluteMomentum(), check_pars.absoluteMomentum(),
                    1.e-8);
  BOOST_CHECK_CLOSE(curv_pars.charge(), check_pars.charge(), 1.e-8);
}

BOOST_AUTO_TEST_CASE(propagator_instatiation_test) {
  auto bField = std::make_shared<NullBField>();
  MultiEigenStepperLoop<> multi_stepper(bField);
  [[maybe_unused]] Propagator<MultiEigenStepperLoop<>> propagator(
      multi_stepper);
}
