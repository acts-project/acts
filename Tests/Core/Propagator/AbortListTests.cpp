// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AbortList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

class Surface;

namespace Test {

  // The path limit abort
  using PathLimit = detail::PathLimitReached;

  // The end of world abort
  using EndOfWorld = detail::EndOfWorldReached;

  // the constrained step class
  using cstep = detail::ConstrainedStep;

  /// This is a simple cache struct to mimic the
  /// Navigator state
  struct NavigatorState
  {
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
  struct PropagatorState
  {

    // This is a simple cache struct to mimic the
    // Stepper cache in the propagation
    struct StepperState
    {
      // accummulated path length cache
      double pathAccumulated = 0.;

      // adaptive sep size of the runge-kutta integration
      cstep stepSize = std::numeric_limits<double>::max();
    };

    /// emulate the options template
    struct Options
    {

      /// The path limit
      double pathLimit = std::numeric_limits<double>::max();

      /// The target tolerance
      double targetTolerance = 1. * units::_um;

      /// Debug output
      /// the string where debug messages are stored (optionally)
      bool        debug       = false;
      std::string debugString = "";
      /// buffer & formatting for consistent output
      size_t debugPfxWidth = 30;
      size_t debugMsgWidth = 50;
    };

    /// Give some options
    Options options;

    /// The stepper state
    StepperState stepping;

    /// The navigator state
    NavigatorState navigation;
  };

  /// This is a struct to mimic the stepper
  struct Stepper
  {
  };

  /// This is a simple result struct to mimic the
  /// propagator result
  struct Result
  {
  };

  // This tests the implementation of the AbortList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(AbortListTest_PathLimit)
  {
    PropagatorState state;
    state.options.pathLimit = 1. * units::_m;
    Stepper stepper;
    Result  result;

    AbortList<PathLimit> abortList;

    // It should not abort yet
    BOOST_CHECK(!abortList(result, state, stepper));
    // The step size should be adapted to 1 meter now
    BOOST_CHECK_EQUAL(state.stepping.stepSize, 1. * units::_m);
    // Let's do a step of 90 cm now
    state.stepping.pathAccumulated = 90. * units::_cm;
    // Still no abort yet
    BOOST_CHECK(!abortList(result, state, stepper));
    // 10 cm are left
    // The step size should be adapted to 10 cm now
    BOOST_CHECK_EQUAL(state.stepping.stepSize, 10. * units::_cm);

    // Approach the target
    while (!abortList(result, state, stepper)) {
      state.stepping.pathAccumulated += 0.5 * state.stepping.stepSize;
    }

    // now we need to be smaller than the tolerance
    BOOST_CHECK_LT(state.stepping.stepSize, 1. * units::_um);

    // Check if you can expand the AbortList
    EndOfWorld eow;
    AbortList<PathLimit, EndOfWorld> pathWorld = abortList.append(eow);
    auto& path = pathWorld.get<PathLimit>();
    BOOST_CHECK(path(result, state, stepper));
  }

}  // namespace Test
}  // namespace Acts
