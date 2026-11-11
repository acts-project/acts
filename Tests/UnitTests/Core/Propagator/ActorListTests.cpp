// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

using namespace Acts;

template <int N>
struct ActorResult {};

template <int N>
struct FailingActor {
  bool fail{};
  bool* called = nullptr;
  using state_type = ActorResult<N>;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                   const navigator_t& /*navigator*/,
                   state_type& /*result*/) const {
    *called = true;
    if (fail) {
      return std::make_error_code(std::errc{});
    } else {
      return Result<void>::success();
    }
  }
};

void testResult(std::array<bool, 3> failing, std::array<bool, 3> expectedCalled,
                bool expectedSuccess) {
  StraightLineStepper stepper{};
  VoidNavigator navigator{};
  Propagator propagator{stepper, navigator};

  std::array<bool, 3> called{false, false, false};
  using ActionList =
      ActorList<FailingActor<1>, FailingActor<2>, FailingActor<3>>;
  Propagator<StraightLineStepper>::Options<ActionList> options{
      Acts::GeometryContext::dangerouslyDefaultConstruct(),
      Acts::MagneticFieldContext{}};

  options.actor.get<FailingActor<1>>().fail = failing[0];
  options.actor.get<FailingActor<1>>().called = &called[0];

  options.actor.get<FailingActor<2>>().fail = failing[1];
  options.actor.get<FailingActor<2>>().called = &called[1];

  options.actor.get<FailingActor<3>>().fail = failing[2];
  options.actor.get<FailingActor<3>>().called = &called[2];

  auto state = propagator.makeState(options);

  Result<void> res =
      options.actor.act(state, stepper, navigator, state.actorState);
  BOOST_CHECK_EQUAL(res.ok(), expectedSuccess);
  BOOST_CHECK_EQUAL(called[0], expectedCalled[0]);
  BOOST_CHECK_EQUAL(called[1], expectedCalled[1]);
  BOOST_CHECK_EQUAL(called[2], expectedCalled[2]);
}

BOOST_AUTO_TEST_CASE(test_all_actors_success) {
  testResult({false, false, false}, {true, true, true}, true);
}

BOOST_AUTO_TEST_CASE(test_first_actor_fails) {
  testResult({true, false, false}, {true, false, false}, false);
}

BOOST_AUTO_TEST_CASE(test_second_actor_fails) {
  testResult({false, true, false}, {true, true, false}, false);
}

BOOST_AUTO_TEST_CASE(test_third_actor_fails) {
  testResult({false, false, true}, {true, true, true}, false);
}

BOOST_AUTO_TEST_CASE(test_first_and_second_actors_fail) {
  testResult({true, true, false}, {true, false, false}, false);
}

BOOST_AUTO_TEST_CASE(all_actors_fail) {
  testResult({true, true, true}, {true, false, false}, false);
}

template <int N>
struct AbortingActor {
  bool abort{};
  bool* called = nullptr;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                  const navigator_t& /*navigator*/) const {
    *called = true;
    return abort;
  }
};

// The composite aborts if any actor aborts, and stops asking at the first one
BOOST_AUTO_TEST_CASE(composite_abort) {
  StraightLineStepper stepper{};
  VoidNavigator navigator{};
  Propagator propagator{stepper, navigator};

  using Actors = ActorList<AbortingActor<1>, AbortingActor<2>>;
  Propagator<StraightLineStepper>::Options<Actors> options{
      Acts::GeometryContext::dangerouslyDefaultConstruct(),
      Acts::MagneticFieldContext{}};

  std::array<bool, 2> called{false, false};
  options.actor.get<AbortingActor<1>>().called = &called[0];
  options.actor.get<AbortingActor<2>>().called = &called[1];

  auto state = propagator.makeState(options);
  BOOST_CHECK(
      !options.actor.checkAbort(state, stepper, navigator, state.actorState));
  BOOST_CHECK(called[0] && called[1]);

  called = {false, false};
  options.actor.get<AbortingActor<1>>().abort = true;
  BOOST_CHECK(
      options.actor.checkAbort(state, stepper, navigator, state.actorState));
  BOOST_CHECK(called[0] && !called[1]);
}
