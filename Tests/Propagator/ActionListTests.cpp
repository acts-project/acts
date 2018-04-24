// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE ActionList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/standard_abort_conditions.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using namespace propagation;

namespace Test {

  // This is a simple cache struct to mimic the
  // Stepper cache in the propagation
  struct Cache
  {
    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive sep size of the runge-kutta integration
    double step_size = std::numeric_limits<double>::max();
  };

  struct DistanceObserver
  {
    double path_to_go = 100. * units::_mm;

    struct this_result
    {
      double distance = std::numeric_limits<double>::max();
    };

    typedef this_result result_type;

    DistanceObserver(double ptg = 0.) : path_to_go(ptg) {}

    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      result.distance = path_to_go - cache.accumulated_path;
    }

    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
    }
  };

  struct CallCounter
  {

    struct this_result
    {
      size_t calls = 0;
    };

    typedef this_result result_type;

    CallCounter() {}

    template <typename cache_t>
    void
    operator()(cache_t&, result_type& r) const
    {
      ++r.calls;
    }

    template <typename cache_t>
    void
    operator()(cache_t&) const
    {
    }
  };

  // This tests the implementation of the ActionList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(ActionListTest_Distance)
  {
    // construct the cache and result
    Cache cache;

    // Type of track parameters produced at the end of the propagation
    typedef typename DistanceObserver::result_type distance_result;
    detail::Extendable<distance_result>            result;

    ActionList<DistanceObserver> action_list;
    action_list.get<DistanceObserver>().path_to_go = 100. * units::_mm;

    // observe and check
    action_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance,
                      100. * units::_mm);

    // now move the cache and check again
    cache.accumulated_path = 50. * units::_mm;
    action_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 50. * units::_mm);
  }

  // This tests the implementation of the ActionList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(ActionListTest_TwoActions)
  {
    // construct the cache and result
    Cache cache;

    // Type of track parameters produced at the end of the propagation
    typedef typename DistanceObserver::result_type distance_result;
    typedef typename CallCounter::result_type      caller_result;
    ActionList<DistanceObserver, CallCounter> action_list;
    action_list.get<DistanceObserver>().path_to_go = 100. * units::_mm;

    detail::Extendable<distance_result, caller_result> result;

    //// observe and check
    action_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance,
                      100. * units::_mm);
    BOOST_CHECK_EQUAL(result.get<caller_result>().calls, 1);

    // now move the cache and check again
    cache.accumulated_path = 50. * units::_mm;
    action_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 50. * units::_mm);
    BOOST_CHECK_EQUAL(result.get<caller_result>().calls, 2);
  }

}  // namespace Test
}  // namespace Acts
