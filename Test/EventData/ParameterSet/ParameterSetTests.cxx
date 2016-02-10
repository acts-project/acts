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

BOOST_AUTO_TEST_CASE(Consistency_checks)
{
  BOOST_CHECK((ParSet_t<ParDefs::loc1,ParDefs::loc2,ParDefs::theta>::size() == 3));
}
/// @endcond
