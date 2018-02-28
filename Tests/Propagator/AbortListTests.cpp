// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
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

#include "ACTS/Propagator/AbortList.hpp"
#include "ACTS/Propagator/detail/Extendable.hpp"
#include "ACTS/Propagator/detail/standard_abort_conditions.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // The path limit abort
  typedef detail::path_limit_reached path_limit;

  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct Cache
  {

    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive sep size of the runge-kutta integration
    double step_size = std::numeric_limits<double>::max();
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