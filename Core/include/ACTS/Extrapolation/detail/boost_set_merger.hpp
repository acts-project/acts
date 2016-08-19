#ifndef ACTS_BOOST_SET_MERGER_HPP
#define ACTS_BOOST_SET_MERGER_HPP 1

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
  struct to_boost_set;

  template <typename... args>
  struct to_boost_set<tuple<args...>>
  {
    typedef typename set<args...>::type type;
  };

  template <typename T, typename R>
  struct to_std_tuple;

  template <typename... args, typename next>
  struct to_std_tuple<tuple<args...>, next>
  {
    typedef tuple<next, args...> type;
  };

  template <typename s, typename v>
  struct boost_set_merger
  {
    typedef typename fold<v, s, insert<_1, _2>>::type unique_types;
    typedef typename fold<unique_types, tuple<>, to_std_tuple<_1, _2>>::type
                                                 flatten;
    typedef typename to_boost_set<flatten>::type type;
  };

}  // namespace detail

}  // namespace Acts

#endif  // ACTS_BOOST_SET_MERGER_HPP
