#ifndef ACTS_ABORT_LIST_CREATOR_HPP
#define ACTS_ABORT_LIST_CREATOR_HPP 1

#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/detail/boost_mpl_helper.hpp"
#include "ACTS/Extrapolation/detail/type_collector.hpp"

namespace Acts {

template <template <typename...> class R, typename... conditions>
struct AbortList : public conditions...
{
  typedef typename detail::type_collector<detail::result_type_extractor,
                                          conditions...>::type abort_results;

  typedef typename detail::unpack_boost_set_as_tparams<R, abort_results>::type
      result_type;
};

}  // namespace Acts

#endif  // ACTS_ABORT_LIST_CREATOR_HPP
