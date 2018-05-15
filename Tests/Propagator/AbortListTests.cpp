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
#include "Acts/Propagator/detail/Extendable.hpp"
#include "Acts/Propagator/detail/standard_abort_conditions.hpp"
#include "Acts/Propagator/detail/constrained_step.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

class Surface;

namespace Test {

  // The path limit abort
  typedef detail::PathLimitReached path_limit;
  // the constrained step class
  typedef detail::ConstrainedStep cstep;

  /// This is a simple cache struct to mimic the
  /// Propagator cache
  struct PropagatorState
  {

    // This is a simple cache struct to mimic the
    // Stepper cache in the propagation
    struct StepperState
    {
      // accummulated path length cache
      double accumulatedPath = 0.;

      // adaptive sep size of the runge-kutta integration
      cstep stepSize = std::numeric_limits<double>::max();
    };

    /// emulate the options template
    struct Options
    {
      /// Debug output
      /// the string where debug messages are stored (optionally)
      bool        debug       = false;
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
    bool           targetReached = false;

    /// Give some options
    Options options;

    /// The Stepper cache
    StepperState stepping;
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
    Result          result;

    AbortList<path_limit> abort_list;

    // Configure path limit with tolerance
    auto& limit           = abort_list.get<path_limit>();
    limit.signedPathLimit = 1. * units::_m;
    limit.tolerance       = 1. * units::_um;

    // It should not abort yet
    BOOST_CHECK(!abort_list(result, state));
    // The step size should be adapted to 1 meter now
    BOOST_CHECK_EQUAL(state.stepping.stepSize, 1. * units::_m);
    // Let's do a step of 90 cm now
    state.stepping.accumulatedPath = 90. * units::_cm;
    // Still no abort yet
    BOOST_CHECK(!abort_list(result, state));
    // 10 cm are left
    // The step size should be adapted to 10 cm now
    BOOST_CHECK_EQUAL(state.stepping.stepSize, 10. * units::_cm);

    // approach the
    while (!abort_list(result, state)) {
      state.stepping.accumulatedPath += 0.5 * state.stepping.stepSize;
    }

    // now we need to be smaller than the tolerance
    BOOST_CHECK(state.stepping.stepSize < 1. * units::_um);
  }

}  // namespace Test
}  // namespace Acts
