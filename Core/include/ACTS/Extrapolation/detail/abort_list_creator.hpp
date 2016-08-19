#ifndef ACTS_ABORT_LIST_CREATOR_HPP
#define ACTS_ABORT_LIST_CREATOR_HPP 1

#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/detail/boost_mpl_helper.hpp"
#include "ACTS/Extrapolation/detail/type_collector.hpp"

namespace Acts {

namespace detail {

  template <template <typename...> class R, int bitmask>
  struct AbortList;

  template <template <typename...> class R, int allbits, int bitmask>
  using AC = AbortCondition<AbortList<R, allbits>, bitmask>;

  template <typename R, int bitmask, typename... bases>
  struct abort_list_creator_impl
  {
    struct AbortListImpl : public bases...
    {
      typedef R result_type;
    };

    typedef AbortListImpl type;
  };

  template <typename R, int bitmask, typename base_sequence>
  struct abort_list_creator_impl_caller;

  template <typename R, int bitmask, typename... bases>
  struct abort_list_creator_impl_caller<R, bitmask, std::tuple<bases...>>
  {
    typedef typename abort_list_creator_impl<R, bitmask, bases...>::type type;
  };

  template <template <typename...> class R, int bitmask>
  struct abort_list_creator
  {
    template <int mask>
    using tAC = AC<R, bitmask, mask>;

    typedef typename type_collector<tAC, bitmask, trait_type_extractor>::type
        active_abort_traits;
    typedef typename type_collector<tAC, bitmask, result_type_extractor>::type
        active_abort_results;

    typedef typename unpack_boost_set_as_tparams<R, active_abort_results>::type
        result_type;

    typedef
        typename abort_list_creator_impl_caller<result_type,
                                                bitmask,
                                                typename boost_set2tuple<active_abort_traits>::
                                                    type>::type abort_list_type;
  };

  template <template <typename...> class R, int bitmask>
  struct AbortList : public abort_list_creator<R, bitmask>::abort_list_type
  {
  };
}  // namespace detail

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_CREATOR_HPP
