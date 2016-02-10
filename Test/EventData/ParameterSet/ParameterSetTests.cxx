/**
 * @file ParameterSetTests.cxx
 *
 * This is a test file.
 */

// Boost include(s)
#define BOOST_TEST_MODULE ParameterSet Tests
#include <boost/test/included/unit_test.hpp>

// ATS include(s)
#include "ParameterSet/ParameterSet.h"
using namespace Ats;

// test include(s)
#include "../../../Test/EventData/ParameterSet/DummyParameterDefs.h"

template<typename ParPolicy::par_id_type... params>
using ParSet_t = ParameterSet<ParPolicy,params...>;

/**
 * @brief Unit test for Ats::anonymous_namespace{ParameterSet.h}::are_sorted helper
 *
 * The test checks for correct behavior in the following cases (always using @c int
 * as value type):
 * -# test: ordered strictly ascending, input: ordered strictly ascending
 * -# test: ordered strictly ascending, input: unordered
 * -# test: ordered strictly ascending, input: ordered weakly ascending
 * -# test: ordered weakly ascending, input: ordered strictly ascending
 * -# test: ordered weakly ascending, input: unordered
 * -# test: ordered weakly ascending, input: ordered weakly ascending
 * -# test: ordered strictly descending, input: ordered strictly descending
 * -# test: ordered strictly descending, input: unordered
 * -# test: ordered strictly descending, input: ordered weakly descending
 * -# test: ordered weakly descending, input: ordered strictly descending
 * -# test: ordered weakly descending, input: unordered
 * -# test: ordered weakly descending, input: ordered weakly descending
 */
BOOST_AUTO_TEST_CASE(are_sorted_helper_tests)
{
  // strictly ascending
  BOOST_CHECK((are_sorted<true,true,int,-1,3,4,12>::value));
  BOOST_CHECK((not are_sorted<true,true,int,-1,13,4>::value));
  BOOST_CHECK((not are_sorted<true,true,int,-1,4,4,7>::value));
  // weakly ascending
  BOOST_CHECK((are_sorted<true,false,int,-1,3,4,12>::value));
  BOOST_CHECK((not are_sorted<true,false,int,-1,13,4>::value));
  BOOST_CHECK((are_sorted<true,false,int,-1,4,4,7>::value));
  // strictly descending
  BOOST_CHECK((are_sorted<false,true,int,1,-3,-4,-12>::value));
  BOOST_CHECK((not are_sorted<false,true,int,1,-13,-4>::value));
  BOOST_CHECK((not are_sorted<false,true,int,1,-4,-4>::value));
  // weakly descending
  BOOST_CHECK((are_sorted<false,false,int,1,-3,-4,-12>::value));
  BOOST_CHECK((not are_sorted<false,false,int,-1,-13,-4>::value));
  BOOST_CHECK((are_sorted<false,false,int,-1,-4,-4,-7>::value));
}

/**
 * @brief Unit test for Ats::anonymous_namespace{ParameterSet.h}::are_within helper
 *
 * The test checks for correct behavior in the following cases (always using @c int
 * as value type):
 * -# all values within (MIN,MAX)
 * -# all values within [MIN,MAX)
 * -# one value < MIN
 * -# multiple values < MIN
 * -# one value > MAX
 * -# multiple values > Max
 * -# one value == MAX
 * -# contains values < MIN and >= MAX
 */
BOOST_AUTO_TEST_CASE(are_within_helper_tests)
{
  BOOST_CHECK((are_within<int,0,10,1,3,7,2>::value));
  BOOST_CHECK((are_within<int,0,10,1,3,0,2>::value));
  BOOST_CHECK((not are_within<int,0,10,-1,3,7,2>::value));
  BOOST_CHECK((not are_within<int,0,10,-1,3,7,-2>::value));
  BOOST_CHECK((not are_within<int,0,10,1,3,17,2>::value));
  BOOST_CHECK((not are_within<int,0,10,1,3,17,12>::value));
  BOOST_CHECK((not are_within<int,0,10,1,10>::value));
  BOOST_CHECK((not are_within<int,0,10,1,-2,10,14>::value));
}
/// @endcond
