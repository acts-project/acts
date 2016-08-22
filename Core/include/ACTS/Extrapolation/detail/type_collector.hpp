#ifndef ACTS_TYPE_COLLECTOR_HPP
#define ACTS_TYPE_COLLECTOR_HPP 1

#include <boost/mpl/set.hpp>
#include "ACTS/Extrapolation/detail/boost_mpl_helper.hpp"
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
    template <template <typename> class extractor,
              typename sequence,
              typename... traits>
    struct impl;

    template <template <typename> class extractor,
              typename sequence,
              typename first,
              typename... others>
    struct impl<extractor, sequence, first, others...>
    {
      typedef typename extractor<first>::type new_type;
      typedef typename impl<extractor,
                            typename bm::insert<sequence, new_type>::type,
                            others...>::type type;
    };

    template <template <typename> class extractor,
              typename sequence,
              typename last>
    struct impl<extractor, sequence, last>
    {
      typedef typename extractor<last>::type new_type;
      typedef typename bm::insert<sequence, new_type>::type type;
    };
  }

  template <template <typename> class extractor, typename... traits>
  struct type_collector
  {
    typedef typename impl<extractor, bm::set<>, traits...>::type found;
    typedef typename boost_set_merger<found, bm::set<>>::type type;
  };
}  // namespace detail

}  // namespace Acts

#endif  // ACTS_TYPE_COLLECTOR_HPP
