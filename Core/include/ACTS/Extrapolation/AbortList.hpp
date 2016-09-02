#ifndef ACTS_ABORT_LIST_CREATOR_HPP
#define ACTS_ABORT_LIST_CREATOR_HPP 1

#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

template <template <typename...> class R, typename... conditions>
struct AbortList : public conditions...
{
private:
  typedef detail::type_collector_t<detail::result_type_extractor, conditions...>
      results;
  typedef detail::type_collector_t<detail::observer_type_extractor,
                                   conditions...>
      observers;

  template <typename... observers>
  using RObserverList = ObserverList<R, observers...>;

public:
  typedef detail::unpack_boost_set_as_tparams_t<R, results> result_type;
  typedef detail::unpack_boost_set_as_tparams_t<RObserverList, observers>
      observer_list_type;
};

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_CREATOR_HPP
