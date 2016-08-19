#ifndef ACTS_BOOST_MPL_HELPER_HPP
#define ACTS_BOOST_MPL_HELPER_HPP 1

#include <boost/mpl/fold.hpp>
#include <boost/mpl/inserter.hpp>
#include <boost/mpl/set.hpp>
#include <tuple>

namespace Acts {

namespace detail {

  namespace {
    using std::tuple;
    using boost::mpl::set;
    using boost::mpl::fold;
    using boost::mpl::insert;
    using boost::mpl::inserter;
    using boost::mpl::placeholders::_1;
    using boost::mpl::placeholders::_2;
  }

  template <typename seq>
  struct tuple2boost_set;

  template <typename... args>
  struct tuple2boost_set<tuple<args...>>
  {
    typedef typename set<args...>::type type;
  };

  template <typename T, typename R>
  struct fold2tuple;

  template <typename... args, typename next>
  struct fold2tuple<tuple<args...>, next>
  {
    typedef tuple<next, args...> type;
  };

  template <typename S>
  struct boost_set2tuple
  {
    typedef typename fold<S, tuple<>, fold2tuple<_1, _2>>::type type;
  };

  template <typename s, typename v>
  struct boost_set_merger
  {
    typedef typename fold<v, s, insert<_1, _2>>::type unique_types;
    typedef
        typename fold<unique_types, tuple<>, fold2tuple<_1, _2>>::type flatten;
    typedef typename tuple2boost_set<flatten>::type type;
  };

  template <template <typename...> class R, typename S>
  struct unpack_boost_set_as_tparams
  {
    typedef typename boost_set2tuple<S>::type as_tuple;
    typedef typename unpack_boost_set_as_tparams<R, as_tuple>::type type;
  };

  template <template <typename...> class R, typename... args>
  struct unpack_boost_set_as_tparams<R, tuple<args...>>
  {
    typedef R<args...> type;
  };

}  // namespace detail

}  // namespace Acts

#endif  // ACTS_BOOST_MPL_HELPER_HPP
