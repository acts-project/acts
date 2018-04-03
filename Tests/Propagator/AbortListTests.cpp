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

namespace Test {

  // The path limit abort
  typedef detail::path_limit_reached path_limit;
  typedef detail::constrained_step   cstep;

  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct Cache
  {

    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive sep size of the runge-kutta integration
    cstep step_size = std::numeric_limits<double>::max();

    /// Debug output
    /// the string where things are stored (optionally)
    std::string debug_string = "";
    /// buffer & formatting for consistent output
    size_t debug_pfx_width = 30;
    size_t debug_msg_width = 50;
    /// flush indication set by actors
    bool debug_flush = false;
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
    Cache  cache;
    Result result;

    AbortList<path_limit> abort_list;

    // Configure path limit with tolerance
    auto& limit             = abort_list.get<path_limit>();
    limit.signed_path_limit = 1. * units::_m;
    limit.tolerance         = 1. * units::_um;

    // It should not abort yet
    BOOST_CHECK(!abort_list(result, cache));
    // The step size should be adapted to 1 meter now
    BOOST_CHECK_EQUAL(cache.step_size, 1. * units::_m);
    // Let's do a step of 90 cm now
    cache.accumulated_path = 90. * units::_cm;
    // Still no abort yet
    BOOST_CHECK(!abort_list(result, cache));
    // 10 cm are left
    // The step size should be adapted to 1 meter now
    BOOST_CHECK_EQUAL(cache.step_size, 10. * units::_cm);

    // approach the
    while (!abort_list(result, cache)) {
      cache.accumulated_path += 0.5 * cache.step_size;
    }

    // now we need to be smaller than the tolerance
    BOOST_CHECK(cache.step_size < 1. * units::_um);
  }

}  // namespace Test
}  // namespace Acts
