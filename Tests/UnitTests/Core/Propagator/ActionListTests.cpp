// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// the constrained step class

/// This is a simple cache struct to mimic the
/// Propagator cache
struct PropagatorState {
  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct StepperState {
    // accummulated path length cache
    double pathAccumulated = 0.;

    // adaptive sep size of the runge-kutta integration
    ConstrainedStep stepSize;
  };

  /// emulate the options template
  struct Options {
    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;
  };

  /// Navigation cache: the start surface
  const Surface* startSurface = nullptr;

  /// Navigation cache: the current surface
  const Surface* currentSurface = nullptr;

  /// Navigation cache: the target surface
  const Surface* targetSurface = nullptr;
  bool targetReached = false;

  /// Give some options
  Options options;

  /// The Stepper cache
  StepperState stepping;
};

/// This is a simple struct to mimic the stepper
struct Stepper {};

/// A distance observer struct as an actor
struct DistanceObserver {
  double path_to_go = 100_mm;

  struct this_result {
    double distance = std::numeric_limits<double>::max();
  };

  using result_type = this_result;

  DistanceObserver(double ptg = 0.) : path_to_go(ptg) {}

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& /*unused*/,
                  result_type& result) const {
    result.distance = path_to_go - state.stepping.pathAccumulated;
  }

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/,
                  const stepper_t& /*unused*/) const {}
};

/// A call counter struct as an actor
struct CallCounter {
  struct this_result {
    size_t calls = 0;
  };

  using result_type = this_result;

  CallCounter() = default;

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/, const stepper_t& /*unused*/,
                  result_type& r) const {
    ++r.calls;
  }

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/,
                  const stepper_t& /*unused*/) const {}
};

// This tests the implementation of the ActionList
// and the standard aborters
BOOST_AUTO_TEST_CASE(ActionListTest_Distance) {
  // construct the (prop/step) cache and result
  PropagatorState state;
  Stepper stepper;

  // Type of track parameters produced at the end of the propagation
  using distance_result = typename DistanceObserver::result_type;
  detail::Extendable<distance_result> result;

  ActionList<DistanceObserver> action_list;
  action_list.get<DistanceObserver>().path_to_go = 100_mm;

  // observe and check
  action_list(state, stepper, result);
  BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 100_mm);

  // now move the cache and check again
  state.stepping.pathAccumulated = 50_mm;
  action_list(state, stepper, result);
  BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 50_mm);
}

// This tests the implementation of the ActionList
// and the standard aborters
BOOST_AUTO_TEST_CASE(ActionListTest_TwoActions) {
  // construct the (prop/step) cache and result
  PropagatorState state;
  Stepper stepper;

  // Type of track parameters produced at the end of the propagation
  using distance_result = typename DistanceObserver::result_type;
  using caller_result = typename CallCounter::result_type;
  ActionList<DistanceObserver, CallCounter> action_list;
  action_list.get<DistanceObserver>().path_to_go = 100_mm;

  detail::Extendable<distance_result, caller_result> result;

  //// observe and check
  action_list(state, stepper, result);
  BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 100_mm);
  BOOST_CHECK_EQUAL(result.get<caller_result>().calls, 1u);

  // now move the cache and check again
  state.stepping.pathAccumulated = 50_mm;
  action_list(state, stepper, result);
  BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 50_mm);
  BOOST_CHECK_EQUAL(result.get<caller_result>().calls, 2u);
}

}  // namespace Test
}  // namespace Acts
