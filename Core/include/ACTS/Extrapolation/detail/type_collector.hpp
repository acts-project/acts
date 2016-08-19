#ifndef ACTS_TYPE_COLLECTOR_HPP
#define ACTS_TYPE_COLLECTOR_HPP 1

#include <boost/mpl/set.hpp>
#include "ACTS/Extrapolation/detail/boost_set_merger.hpp"
namespace bm = boost::mpl;

namespace Acts {

namespace detail {

  template <typename T>
  struct result_type_extractor
  {
    typedef typename T::result_type type;
  };

  template <typename T>
  struct observer_type_extractor
  {
    typedef typename T::observer_type type;
  };

  template <typename T>
  struct trait_type_extractor
  {
    typedef T type;
  };

  namespace {
    template <template <int> class traits,
              int bitmask,
              template <typename> class extractor,
              typename sequence>
    struct impl;

    template <template <int> class traits,
              int bitmask,
              template <typename> class extractor,
              typename sequence>
    struct impl
    {
      static constexpr int lsb        = bitmask & ~(bitmask - 1);
      static constexpr int other_bits = bitmask ^ lsb;

      typedef typename extractor<traits<lsb>>::type new_type;
      typedef typename impl<traits,
                            other_bits,
                            extractor,
                            typename bm::insert<sequence, new_type>::type>::type
          type;
    };

    template <template <int> class traits,
              template <typename> class extractor,
              typename sequence>
    struct impl<traits, 0, extractor, sequence>
    {
      typedef sequence type;
    };
  }

  template <template <int> class traits,
            int bitmask,
            template <typename> class extractor>
  struct type_collector
  {
    typedef typename impl<traits, bitmask, extractor, bm::set<>>::type found;
    typedef typename boost_set_merger<found, bm::set<>>::type type;
  };
}  // namespace detail

}  // namespace Acts

#endif  // ACTS_TYPE_COLLECTOR_HPP
