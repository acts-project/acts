// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {

class Surface;

namespace Test {

// The path limit abort
using PathLimit = PathLimitReached;

// The end of world abort
using EndOfWorld = EndOfWorldReached;

// the constrained step class

/// This is a simple cache struct to mimic the
/// Navigator state
struct NavigatorState {
  /// Navigation cache: the start surface
  const Surface* startSurface = nullptr;

  /// Navigation cache: the current surface
  const Surface* currentSurface = nullptr;

  /// Navigation cache: the target surface
  const Surface* targetSurface = nullptr;

  /// The boolean is reached
  bool targetReached = false;
};

/// This is a simple cache struct to mimic the
/// Propagator state
struct PropagatorState {
  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct StepperState {
    // Accummulated path length cache
    double pathAccumulated = 0.;
    // Navigation direction
    NavigationDirection navDir = NavigationDirection::Forward;
    // adaptive sep size of the runge-kutta integration
    ConstrainedStep stepSize;
  };

  /// emulate the options template
  struct Options {
    /// The path limit
    double pathLimit = std::numeric_limits<double>::max();

    /// The target tolerance
    double targetTolerance = 1_um;

    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = true;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;

    LoggerWrapper logger{getDummyLogger()};
  };

  /// Give some options
  Options options;

  /// The stepper state
  StepperState stepping;

  /// The navigator state
  NavigatorState navigation;
};

/// This is a struct to mimic the stepper
struct Stepper {
  auto outputStepSize(const PropagatorState::StepperState&) const {
    return std::string{};
  }

  double getStepSize(const PropagatorState::StepperState& state,
                     ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  void setStepSize(PropagatorState::StepperState& state, double stepSize,
                   ConstrainedStep::Type stype = ConstrainedStep::actor,
                   bool release = true) const {
    state.stepSize.update(stepSize, stype, release);
  }
};

/// This is a simple result struct to mimic the
/// propagator result
struct Result {};

// This tests the implementation of the AbortList
// and the standard aborters
BOOST_AUTO_TEST_CASE(AbortListTest_PathLimit) {
  PropagatorState state;
  state.options.pathLimit = 1_m;
  Stepper stepper;
  Result result;

  AbortList<PathLimit> abortList;
  abortList.get<PathLimit>().internalLimit = state.options.pathLimit;

  // It should not abort yet
  BOOST_CHECK(!abortList(result, state, stepper));
  // The step size should be adapted to 1 meter now
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), 1_m);
  // Let's do a step of 90 cm now
  state.stepping.pathAccumulated = 90_cm;
  // Still no abort yet
  BOOST_CHECK(!abortList(result, state, stepper));
  // 10 cm are left
  // The step size should be adapted to 10 cm now
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), 10_cm);

  // Approach the target
  while (!abortList(result, state, stepper)) {
    state.stepping.pathAccumulated += 0.5 * state.stepping.stepSize.value();
  }

  // now we need to be smaller than the tolerance
  BOOST_CHECK_LT(state.stepping.stepSize.value(), 1_um);

  // Check if you can expand the AbortList
  EndOfWorld eow;
  AbortList<PathLimit, EndOfWorld> pathWorld = abortList.append(eow);
  auto& path = pathWorld.get<PathLimit>();
  BOOST_CHECK(path(state, stepper));
}

struct ActorA {
  struct result_type {
    int a_number{42};
  };

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t&, const stepper_t&, result_type&) const {}
};

struct ActorB {
  struct result_type {
    int a_number{42};
  };

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t&, const stepper_t&, result_type&) const {}
};

struct AborterWithResultA {
  using action_type = ActorA;

  template <typename propagator_state_t, typename stepper_t, typename result_t>
  bool operator()(const propagator_state_t&, const stepper_t&,
                  const result_t&) const {
    return true;
  }
};

struct AborterWithResultInvalid {
  using action_type = ActorA;

  template <typename propagator_state_t, typename stepper_t>
  bool operator()(const propagator_state_t&, const stepper_t&) const {
    return true;
  }
};

struct AborterWithResultB {
  using action_type = ActorB;

  template <typename propagator_state_t, typename stepper_t, typename result_t>
  bool operator()(const propagator_state_t&, const stepper_t&,
                  const result_t&) const {
    return true;
  }

  template <typename propagator_state_t, typename stepper_t>
  bool operator()(const propagator_state_t&, const stepper_t&) const {
    return true;
  }
};

struct AborterWithoutResult {
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(const propagator_state_t&, const stepper_t&) const {
    return true;
  }
};

struct AborterWithoutResultInvalid {
  template <typename propagator_state_t, typename stepper_t, typename result_t>
  bool operator()(const propagator_state_t&, const stepper_t&,
                  const result_t&) const {
    return true;
  }
};

template <typename P, typename S, typename... As>
constexpr bool signature_check =
    detail::all_of_v<Concepts ::abort_condition_signature_check_v<As, P, S>...>;

BOOST_AUTO_TEST_CASE(AbortListSignatureTest) {
  using P = PropagatorState;
  using S = Stepper;

  static_assert(signature_check<P, S, AborterWithoutResult>, "failed");
  static_assert(signature_check<P, S, AborterWithResultA>, "failed");
  static_assert(signature_check<P, S, AborterWithResultB>, "failed");

  // combination of two valid aborters
  static_assert(signature_check<P, S, AborterWithoutResult, AborterWithResultA,
                                AborterWithResultB>,
                "failed");

  // declares an `action_type`, but call operator is invalid since it has no
  // result argument
  static_assert(!signature_check<P, S, AborterWithResultInvalid>, "failed");
  // does not declare an `action_type` but the call operator expects a result
  // argument
  static_assert(!signature_check<P, S, AborterWithoutResultInvalid>, "failed");
}

}  // namespace Test

}  // namespace Acts
