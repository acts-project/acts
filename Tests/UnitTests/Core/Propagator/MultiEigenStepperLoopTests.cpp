// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperDefaultExtension.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"

#include <memory>
#include <optional>

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

BOOST_AUTO_TEST_CASE(
    cylindrical_portal_guard_preserves_multi_stepper_navigation) {
  using namespace Acts::UnitLiterals;
  const auto& gctx = tester.geoCtx;
  const auto& mctx = tester.magCtx;
  Blueprint::Config config;
  config.envelope[AxisDirection::AxisR] = {100_mm, 100_mm};
  config.envelope[AxisDirection::AxisZ] = {100_mm, 100_mm};
  Blueprint root(config);
  root.addStaticVolume(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(100_mm, 200_mm, 500_mm),
      "annulus");
  std::shared_ptr<const TrackingGeometry> geometry = root.construct({}, gctx);
  Propagator<MultiEigenStepperLoop<>, Navigator> propagator(
      MultiEigenStepperLoop<>(
          std::make_shared<ConstantBField>(Vector3(0., 0., 2_T))),
      Navigator({.trackingGeometry = geometry}));

  const auto single = BoundTrackParameters::createCurvilinear(
      Vector4(150_mm, 0., 0., 0.), Vector3::UnitX(), 1. / 10_GeV, std::nullopt,
      ParticleHypothesis::pion());
  MultiComponentBoundTrackParameters parameters(
      single.referenceSurface().getSharedPtr(), false,
      ParticleHypothesis::pion());
  parameters.pushComponent(0.5, single.parameters(), std::nullopt);
  parameters.pushComponent(0.5, single.parameters(), std::nullopt);

  decltype(propagator)::Options<> options(gctx, mctx);
  options.pathLimit = 100_mm;
  options.loopProtection = false;
  options.maxSteps = 100;
  options.stepping.maxStepSize = 20_mm;
  options.navigation.experimentalCylindricalPortalGuard = false;
  const auto disabled = propagator.propagate(parameters, options, false);
  BOOST_REQUIRE(disabled.ok());
  BOOST_CHECK_GT(disabled->steps, 1u);
  BOOST_CHECK_SMALL(disabled->pathLength - options.pathLimit,
                    s_onSurfaceTolerance);

  // The optional single-trajectory guard must compile and remain a no-op for
  // composite steppers, which have no single ConstrainedStep accessor.
  options.navigation.experimentalCylindricalPortalGuard = true;
  const auto enabled = propagator.propagate(parameters, options, false);
  BOOST_REQUIRE(enabled.ok());
  BOOST_CHECK_EQUAL(enabled->steps, disabled->steps);
  BOOST_CHECK_EQUAL(enabled->pathLength, disabled->pathLength);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
