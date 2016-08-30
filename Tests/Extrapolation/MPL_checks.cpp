
// Boost include(s)
#define BOOST_TEST_MODULE MPL Tests
#include <ACTS/Extrapolation/AbortList.hpp>
#include <ACTS/Utilities/detail/MPL/all_of.hpp>
#include <ACTS/Utilities/detail/MPL/any_of.hpp>
#include <ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp>
#include <ACTS/Utilities/detail/MPL/has_duplicates.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <type_traits>
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

#include <iostream>

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(boost_set_merger_test)
  {
    typedef typename boost::mpl::set<float, int, char, bool>::type first;
    typedef typename boost::mpl::vector<long, int, void*>::type second;
    typedef typename detail::boost_set_merger_t<first, second> found;
    typedef typename boost::mpl::set<float, int, char, bool, long, void*>::type
        expected;

    static_assert(std::is_same<found, expected>::value,
                  "merging sequence into boost::mpl::set failed");
  }

  template <typename... args>
  struct variadic_struct
  {
  };

  BOOST_AUTO_TEST_CASE(unpack_boost_set_as_template_test)
  {
    typedef boost::mpl::set<float, int, char>::type boost_set;
    typedef variadic_struct<float, int, char>       expected;
    typedef detail::unpack_boost_set_as_tparams_t<variadic_struct, boost_set>
        found;

    static_assert(std::is_same<found, expected>::value,
                  "using boost::mpl::set for variadic templates failed");
  }

  namespace {
    struct traits1
    {
      typedef int  result_type;
      typedef char observer_type;
    };

    template <bool>
    struct traits2;

    template <>
    struct traits2<false>
    {
      typedef bool  result_type;
      typedef float observer_type;
    };

    template <>
    struct traits2<true>
    {
      typedef float observer_type;
    };
  }

  BOOST_AUTO_TEST_CASE(type_collector_test)
  {
    typedef detail::type_collector_t<detail::result_type_extractor,
                                     traits1,
                                     traits2<true>,
                                     traits2<false>>
        found_results;

    typedef detail::type_collector_t<detail::observer_type_extractor,
                                     traits1,
                                     traits2<true>,
                                     traits2<false>>
        found_observers;

    typedef typename boost::mpl::set<int, bool>::type   expected_results;
    typedef typename boost::mpl::set<char, float>::type expected_observers;

    static_assert(std::is_same<found_results, expected_results>::value,
                  "collecting result types failed");
    static_assert(std::is_same<found_observers, expected_observers>::value,
                  "collecting observer types failed");
  }

  BOOST_AUTO_TEST_CASE(has_duplicates_test)
  {
    using detail::has_duplicates_v;
    static_assert(has_duplicates_v<int, float, char, int>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, int, char, float>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, char, float, float>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, char, char, float>,
                  "has_duplicates_v failed");
    static_assert(not has_duplicates_v<int, bool, char, float>,
                  "has_duplicates_v failed");
  }

  BOOST_AUTO_TEST_CASE(all_of_test)
  {
    using detail::all_of_v;

    static_assert(not all_of_v<true, true, false>,
                  "all_of_v<true, true, false> failed");
    static_assert(not all_of_v<false, true, true, false>,
                  "all_of_v<false, true, true, false> failed");
    static_assert(all_of_v<true, true, true>,
                  "all_of_v<true, true, true> failed");
    static_assert(all_of_v<true>, "all_of_v<true> failed");
    static_assert(not all_of_v<false>, "all_of_v<false> failed");
    static_assert(all_of_v<>, "all_of_v<> failed");
  }

  BOOST_AUTO_TEST_CASE(any_of_test)
  {
    using detail::any_of_v;

    static_assert(any_of_v<true, true, false>,
                  "any_of_v<true, true, false> failed");
    static_assert(any_of_v<false, true, true, false>,
                  "any_of_v<false, true, true, false> failed");
    static_assert(any_of_v<true, true, true>,
                  "any_of_v<true, true, true> failed");
    static_assert(not any_of_v<false, false>, "any_of_v<false, false> failed");
    static_assert(any_of_v<true>, "any_of_v<true> failed");
    static_assert(not any_of_v<false>, "any_of_v<false> failed");
    static_assert(not any_of_v<>, "any_of_v<> failed");
  }

  BOOST_AUTO_TEST_CASE(abort_list_test)
  {
    typedef AbortList<variadic_struct,
                      DestinationSurface,
                      MaxMaterial,
                      MaxPathLength>
                                      AList;
    typedef AList::result_type        found_result_type;
    typedef AList::observer_list_type found_observer_list_type;
    typedef variadic_struct<MaxMaterial::result_type,
                            MaxPathLength::result_type>
        expected_result_type;

    static_assert(std::is_same<found_result_type, expected_result_type>::value,
                  "creating abort condition list failed");

    static_assert(std::is_base_of<DestinationSurface, AList>::value, "");
    static_assert(std::is_base_of<MaxMaterial, AList>::value, "");
    static_assert(std::is_base_of<MaxPathLength, AList>::value, "");

    AList al;
    al.maxMaterial   = 3;
    al.maxPathLength = 17.8;

    AList c = al;
    std::cout << c.maxMaterial << std::endl;
    std::cout << c.maxPathLength << std::endl;
    std::cout << sizeof(c) << std::endl;
    std::cout << sizeof(double) << std::endl;
  }
}  // namespace Test

}  // namespace Acts
