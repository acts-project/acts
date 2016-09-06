#ifndef ACTS_TYPE_COLLECTOR_HPP
#define ACTS_TYPE_COLLECTOR_HPP 1

#include <boost/mpl/set.hpp>
#include <type_traits>
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
namespace bm = boost::mpl;

#define HAS_TYPE(x)                                                            \
  template <typename T, typename = void>                                       \
  struct has_type : public std::false_type                                     \
  {                                                                            \
  };                                                                           \
                                                                               \
  template <typename T>                                                        \
  struct has_type<T,                                                           \
                  typename std::enable_if<(sizeof(typename T::x) > 0),         \
                                          void>::type> : public std::true_type \
  {                                                                            \
  };

namespace Acts {

namespace detail {

  namespace {

    struct result_type_extractor
    {
    private:
      template <typename T>
      struct extractor
      {
        typedef typename T::result_type type;
      };

    public:
      HAS_TYPE(result_type);

      template <typename T>
      using type = typename extractor<T>::type;
    };

    struct observer_type_extractor
    {
    private:
      template <typename T>
      struct extractor
      {
        typedef typename T::observer_type type;
      };

    public:
      HAS_TYPE(observer_type);

      template <typename T>
      using type = typename extractor<T>::type;
    };

    template <typename sequence, typename ex, typename T, bool has_type = false>
    struct type_inserter
    {
      typedef sequence type;
    };

    template <typename sequence, typename ex, typename T>
    struct type_inserter<sequence, ex, T, true>
    {
      typedef typename bm::insert<sequence, typename ex::template type<T>>::type
          type;
    };

    template <typename sequence, typename ex, typename... traits>
    struct type_collector_impl;

    template <typename sequence, typename ex>
    struct type_collector_impl<sequence, ex>
    {
      typedef sequence type;
    };

    template <typename sequence,
              typename ex,
              typename first,
              typename... others>
    struct type_collector_impl<sequence, ex, first, others...>
    {
      typedef typename type_inserter<sequence,
                                     ex,
                                     first,
                                     ex::template has_type<first>::value>::type
          new_seq;
      typedef typename type_collector_impl<new_seq, ex, others...>::type type;
    };

    template <typename sequence, typename ex, typename last>
    struct type_collector_impl<sequence, ex, last>
    {
      typedef
          typename type_inserter<sequence,
                                 ex,
                                 last,
                                 ex::template has_type<last>::value>::type type;
    };

    template <typename extractor, typename... traits>
    struct type_collector
    {
      typedef
          typename type_collector_impl<bm::set<>, extractor, traits...>::type
              found;
      typedef boost_set_merger_t<found, bm::set<>> type;
    };
  }  // end of anonymous namespace

  template <typename T>
  constexpr bool has_result_type_v = result_type_extractor::has_type<T>::value;

  template <typename T>
  using result_type_t = typename result_type_extractor::type<T>;

  template <typename T>
  constexpr bool has_observer_type_v
      = observer_type_extractor::has_type<T>::value;

  template <typename T>
  using observer_type_t = typename observer_type_extractor::type<T>;

  template <typename extractor, typename... traits>
  using type_collector_t = typename type_collector<extractor, traits...>::type;
}  // namespace detail

}  // namespace Acts

#endif  // ACTS_TYPE_COLLECTOR_HPP
