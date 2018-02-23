// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE ObserverList Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/ObserverList.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Propagator/detail/standard_abort_conditions.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

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

    template <typename input_t>
    void
    operator()(const input_t& cache, result_type& result) const
    {
      result.distance = path_to_go - cache.accumulated_path;
    }

    template <typename input_t>
    void
    operator()(const input_t& cache) const
    {
      (void)cache;
    }
  };

  // This tests teh implementation of the ObserverList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(ObserverListTest_Distance)
  {
    // construct the cache and result
    Cache cache;

    // Type of track parameters produced at the end of the propagation
    typedef typename DistanceObserver::result_type distance_result;
    distance_result                                dr;
    detail::Extendable<distance_result>            result(dr);

    DistanceObserver d100(100. * units::_mm);
    auto             observer_list = ObserverList<DistanceObserver>({d100});

    // observe and check
    observer_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance,
                      100. * units::_mm);

    // now move the cache and check again
    cache.accumulated_path = 50. * units::_mm;
    observer_list(cache, result);
    BOOST_CHECK_EQUAL(result.get<distance_result>().distance, 50. * units::_mm);
  }

}  // namespace Test
}  // namespace Acts