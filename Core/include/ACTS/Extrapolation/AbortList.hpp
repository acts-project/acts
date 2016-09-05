#ifndef ACTS_ABORT_LIST_CREATOR_HPP
#define ACTS_ABORT_LIST_CREATOR_HPP 1

#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

template <template <typename...> class R,
          typename input,
          typename... conditions>
struct AbortList
{
private:
  static_assert(not detail::has_duplicates_v<conditions...>,
                "same abort conditions type specified several times");

private:
  typedef detail::type_collector_t<detail::observer_type_extractor,
                                   conditions...>
      observers;

  template <typename... observers>
  using ThisObserverList = ObserverList<R, input, observers...>;

public:
  typedef detail::unpack_boost_set_as_tparams_t<ThisObserverList, observers>
      observer_list_type;
};

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_CREATOR_HPP
