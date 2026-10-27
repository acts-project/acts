// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperDefaultExtension.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include "MultiStepperTests.hpp"

using MyEigenStepper = EigenStepper<EigenStepperDefaultExtension>;

const MultiStepperTester<MyEigenStepper, MultiStepperLoop<MyEigenStepper>>
    tester;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(multi_stepper_config_constructor) {
  tester.test_config_constructor();
}
BOOST_AUTO_TEST_CASE(multi_stepper_state_no_cov) {
  tester.test_multi_stepper_state<false>();
}

BOOST_AUTO_TEST_CASE(test_combined_bound_state) {
  tester.test_combined_bound_state_function();
}

BOOST_AUTO_TEST_CASE(test_surface_status_and_cmpwise_bound_state) {
  tester.test_multi_stepper_surface_status_update();
}

BOOST_AUTO_TEST_CASE(multi_eigen_vs_single_eigen) {
  tester.test_multi_stepper_vs_eigen_stepper();
}

BOOST_AUTO_TEST_CASE(multi_eigen_component_iterable_with_modification) {
  tester.test_components_modifying_accessors();
}

BOOST_AUTO_TEST_CASE(test_single_component_interface) {
  tester.test_single_component_interface_function();
}

BOOST_AUTO_TEST_CASE(remove_add_components_test) {
  tester.remove_add_components_function();
}

BOOST_AUTO_TEST_CASE(test_component_wise_bound_state) {
  tester.test_component_bound_state();
}

BOOST_AUTO_TEST_CASE(test_curvilinear_state) {
  tester.test_combined_curvilinear_state_function();
}

BOOST_AUTO_TEST_CASE(propagator_instatiation_test) {
  tester.propagator_instatiation_test_function();
}

BOOST_AUTO_TEST_CASE(component_curvature_controls_bounds) {
  using Stepper = MultiStepperLoop<MyEigenStepper>;
  Stepper stepper(std::make_shared<ConstantBField>(Vector3(0.005, 0., 0.)));
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  const auto start =
      CurvilinearSurface(Vector3::Zero(), Vector3::UnitZ()).planeSurface();
  std::vector<std::tuple<double, BoundVector>> components(
      2, {0.5, BoundVector::Zero()});
  std::get<1>(components[0])[eBoundQOverP] = 1.;
  std::get<1>(components[1])[eBoundQOverP] = 0.5;
  MultiComponentBoundTrackParameters parameters(
      start, components,
      [](const auto& component) {
        return std::tie(std::get<0>(component), std::get<1>(component));
      },
      ParticleHypothesis::pion());
  auto state = stepper.makeState(Stepper::Options(gctx, mctx));
  stepper.initialize(state, parameters);
  const auto sensor = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(0., 1., 20.)),
      std::make_shared<RectangleBounds>(0.1, 0.1));
  const auto result = stepper.updateSurfaceStatusWithCurvature(
      state, *sensor, 0, Direction::Forward(), BoundaryTolerance::None(),
      s_onSurfaceTolerance, ConstrainedStep::Type::Navigator);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(*result == IntersectionStatus::reachable);
  BOOST_CHECK(state.components[0].status == IntersectionStatus::reachable);
  BOOST_CHECK(state.components[1].status == IntersectionStatus::unreachable);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
