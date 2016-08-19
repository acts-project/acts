
// Boost include(s)
#define BOOST_TEST_MODULE ParameterSet Tests
#include <ACTS/Extrapolation/detail/boost_mpl_helper.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <type_traits>
#include "ACTS/Extrapolation/detail/abort_list_creator.hpp"
#include "ACTS/Extrapolation/detail/type_collector.hpp"

#include <iostream>

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(boost_set_merger_test)
  {
    typedef typename boost::mpl::set<float, int, char, bool>::type first;
    typedef typename boost::mpl::vector<long, int, void*>::type second;
    typedef typename detail::boost_set_merger<first, second>::type found;
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
    typedef typename detail::unpack_boost_set_as_tparams<variadic_struct,
                                                         boost_set>::type found;

    static_assert(std::is_same<found, expected>::value,
                  "using boost::mpl::set for variadic templates failed");
  }

  namespace {
    template <int bitmask>
    struct dummy_traits;

    template <>
    struct dummy_traits<1 << 0>
    {
      typedef int  result_type;
      typedef char observer_type;
    };

    template <>
    struct dummy_traits<1 << 1>
    {
      typedef bool  result_type;
      typedef float observer_type;
    };

    template <>
    struct dummy_traits<1 << 2>
    {
      typedef struct
      {
      } result_type;
      typedef float observer_type;
    };
  }

  BOOST_AUTO_TEST_CASE(type_collector_test)
  {
    typedef typename detail::type_collector<dummy_traits,
                                            7,
                                            detail::result_type_extractor>::type
        found_results;

    typedef
        typename detail::type_collector<dummy_traits,
                                        7,
                                        detail::observer_type_extractor>::type
            found_observers;

    typedef typename detail::type_collector<dummy_traits,
                                            7,
                                            detail::trait_type_extractor>::type
        found_traits;

    typedef
        typename boost::mpl::set<int, bool, dummy_traits<4>::result_type>::type
            expected_results;
    typedef typename boost::mpl::set<char, float>::type expected_observers;
    typedef typename boost::mpl::set<dummy_traits<1>,
                                     dummy_traits<2>,
                                     dummy_traits<4>>::type expected_traits;

    static_assert(std::is_same<found_results, expected_results>::value,
                  "collecting result types failed");
    static_assert(std::is_same<found_observers, expected_observers>::value,
                  "collecting observer types failed");
    static_assert(std::is_same<found_traits, expected_traits>::value,
                  "collecting trait types failed");
  }

  template <int allbits, int bitmask>
  using AC
      = AbortCondition<detail::AbortList<variadic_struct, allbits>, bitmask>;

  template <typename... Args>
  void
  f()
  {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }

  BOOST_AUTO_TEST_CASE(abort_list_creator_test)
  {
    constexpr int bitmask
        = Acts::DestinationSurface | Acts::MaxMaterial | Acts::MaxPathLength;
    typedef typename detail::abort_list_creator<variadic_struct,
                                                bitmask>::result_type found;
    typedef variadic_struct<AC<bitmask, Acts::DestinationSurface>::result_type,
                            AC<bitmask, Acts::MaxMaterial>::result_type,
                            AC<bitmask, Acts::MaxPathLength>::result_type>
        expected;

    static_assert(std::is_same<found, expected>::value,
                  "creating abort condition list failed");

    static_assert(
        std::is_base_of<AC<bitmask, Acts::DestinationSurface>,
                        detail::AbortList<variadic_struct, bitmask>>::value,
        "");
    static_assert(
        std::is_base_of<AC<bitmask, Acts::MaxMaterial>,
                        detail::AbortList<variadic_struct, bitmask>>::value,
        "");
    static_assert(
        std::is_base_of<AC<bitmask, Acts::MaxPathLength>,
                        detail::AbortList<variadic_struct, bitmask>>::value,
        "");

    f<detail::abort_list_creator<variadic_struct,
                                 bitmask>::active_abort_traits>();

    detail::AbortList<variadic_struct, bitmask> al;
    al.destinationSurface().maxMaterial(3).maxPathLength(17.8);

    std::cout << al.m_maxMaterial << std::endl;
    std::cout << al.m_maxPathLength << std::endl;
  }
}  // namespace Test

}  // namespace Acts
