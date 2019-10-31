// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE FSM Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include "Acts/Utilities/FiniteStateMachine.hpp"

#include <iostream>

namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

namespace states {
struct Disconnected {};

struct Connecting {};
struct Pinging {};
struct Connected {};
}  // namespace states

namespace events {
struct Connect {};
struct Established {};
struct Timeout {};
struct Ping {};
struct Pong {};
struct Disconnect {};
}  // namespace events

struct fsm : FiniteStateMachine<fsm, states::Disconnected, states::Connecting,
                                states::Pinging, states::Connected> {
  fsm() : fsm_base(states::Disconnected{}){};

  event_return on_event(const states::Disconnected&, const events::Connect&) {
    return states::Connecting{};
  }

  event_return on_event(const states::Connecting&, const events::Established&) {
    return states::Connected{};
  }

  event_return on_event(const states::Connected&, const events::Ping&) {
    std::cout << "ping!" << std::endl;
    setState(states::Pinging{});
    return process_event(events::Pong{});
  }

  event_return on_event(const states::Pinging&, const events::Pong&) {
    std::cout << "pong!" << std::endl;
    return states::Connected{};
  }

  event_return on_event(const states::Connected&, const events::Timeout&) {
    return states::Connecting{};
  }

  event_return on_event(const states::Connected&, const events::Disconnect&) {
    return states::Disconnected{};
  }

  template <typename State, typename Event>
  event_return on_event(const State&, const Event&) const {
    return Terminated{};
  }

  template <typename State, typename... Args>
  void on_enter(const State&, Args&&...) {}

  template <typename State, typename... Args>
  void on_exit(const State&, Args&&...) {}

  template <typename... Args>
  void log(Args&&...) {}
};

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(Transitions) {
  fsm sm{};
  BOOST_CHECK(sm.is(states::Disconnected{}));
  sm.dispatch(events::Connect{});
  BOOST_CHECK(sm.is(states::Connecting{}));
  sm.dispatch(events::Established{});
  BOOST_CHECK(sm.is(states::Connected{}));
  sm.dispatch(events::Ping{});
  sm.dispatch(events::Ping{});
  sm.dispatch(events::Ping{});
  sm.dispatch(events::Ping{});
  BOOST_CHECK(sm.is(states::Connected{}));
  sm.dispatch(events::Timeout{});
  BOOST_CHECK(sm.is(states::Connecting{}));
  sm.dispatch(events::Established{});
  BOOST_CHECK(sm.is(states::Connected{}));
  sm.dispatch(events::Disconnect{});
  BOOST_CHECK(sm.is(states::Disconnected{}));
}

BOOST_AUTO_TEST_CASE(Terminted) {
  fsm sm{};
  BOOST_CHECK(sm.is(states::Disconnected{}));

  sm.dispatch(events::Disconnect{});
  BOOST_CHECK(sm.terminated());
}

struct fsm2
    : FiniteStateMachine<fsm2, states::Disconnected, states::Connected> {
  fsm2() : fsm_base(states::Disconnected{}){};

  event_return on_event(const states::Disconnected&, const events::Connect&,
                        double f) {
    std::cout << "f: " << f << std::endl;
    return states::Connected{};
  }

  event_return on_event(const states::Connected&, const events::Disconnect&) {
    std::cout << "disconnect!" << std::endl;
    return states::Disconnected{};
  }

  template <typename State, typename Event, typename... Args>
  event_return on_event(const State&, const Event&, Args&&...) const {
    return Terminated{};
  }

  template <typename... Args>
  void on_enter(const Terminated&, Args&&...) {
    throw std::runtime_error("FSM terminated!");
  }

  template <typename State, typename... Args>
  void on_enter(const State&, Args&&...) {}

  template <typename State, typename... Args>
  void on_exit(const State&, Args&&...) {}
  template <typename... Args>
  void log(Args&&...) {}
};

BOOST_AUTO_TEST_CASE(Arguments) {
  fsm2 sm{};
  BOOST_CHECK(sm.is(states::Disconnected{}));

  sm.dispatch(events::Connect{}, 42.);
  BOOST_CHECK(sm.is(states::Connected{}));
  sm.dispatch(events::Disconnect{});
  BOOST_CHECK(sm.is(states::Disconnected{}));
  sm.dispatch(events::Connect{}, -1.);

  // call disconnect, but disconnect does not accept this call signature
  BOOST_REQUIRE_THROW(sm.dispatch(events::Disconnect{}, 9), std::runtime_error);
  BOOST_CHECK(sm.terminated());

  // cannot dispatch on terminated (in this specific configuration, in
  // general terminated is just another state).
  BOOST_REQUIRE_THROW(sm.dispatch(events::Connect{}), std::runtime_error);
  // still in terminated
  BOOST_CHECK(sm.terminated());

  // we can reset the state though!
  sm.setState(states::Disconnected{});
  BOOST_CHECK(sm.is(states::Disconnected{}));
  sm.dispatch(events::Connect{}, -1.);
  BOOST_CHECK(sm.is(states::Connected{}));
}

struct S1 {};
struct S2 {};
struct S3 {};

struct E1 {};
struct E2 {};
struct E3 {};

struct fsm3 : FiniteStateMachine<fsm3, S1, S2, S3> {
  bool on_exit_called = false;
  bool on_enter_called = false;
  void reset() {
    on_exit_called = false;
    on_enter_called = false;
  }

  // S1 + E1 = S2
  event_return on_event(const S1&, const E1&) { return S2{}; }

  // S2 + E1 = S2
  // external transition to self
  event_return on_event(const S2&, const E1&) { return S2{}; }

  // S2 + E2
  // internal transition
  event_return on_event(const S2&, const E2&) {
    return std::nullopt;
    // return S2{};
  }

  // S2 + E3 = S3
  // external transition
  event_return on_event(const S2&, const E3&) { return S3{}; }

  // catchers

  template <typename State, typename Event, typename... Args>
  event_return on_event(const State&, const Event&, Args&&...) const {
    return Terminated{};
  }

  template <typename State, typename... Args>
  void on_enter(const State&, Args&&...) {
    on_enter_called = true;
  }

  template <typename State, typename... Args>
  void on_exit(const State&, Args&&...) {
    on_exit_called = true;
  }
  template <typename... Args>
  void log(Args&&...) {}
};

BOOST_AUTO_TEST_CASE(InternalTransitions) {
  fsm3 sm;
  BOOST_CHECK(sm.is(S1{}));

  sm.dispatch(E1{});
  BOOST_CHECK(sm.is(S2{}));
  BOOST_CHECK(sm.on_exit_called);
  BOOST_CHECK(sm.on_enter_called);

  sm.reset();

  sm.dispatch(E1{});
  // still in S2
  BOOST_CHECK(sm.is(S2{}));
  // on_enter / exit should have been called
  BOOST_CHECK(sm.on_exit_called);
  BOOST_CHECK(sm.on_enter_called);
  sm.reset();

  sm.dispatch(E2{});
  // still in S2
  BOOST_CHECK(sm.is(S2{}));
  // on_enter / exit should NOT have been called
  BOOST_CHECK(!sm.on_exit_called);
  BOOST_CHECK(!sm.on_enter_called);
  sm.reset();

  sm.dispatch(E3{});
  BOOST_CHECK(sm.is(S3{}));
  // on_enter / exit should have been called
  BOOST_CHECK(sm.on_exit_called);
  BOOST_CHECK(sm.on_enter_called);
  sm.reset();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
