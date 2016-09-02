#include <boost/mpl/set.hpp>
#include <functional>
#include <vector>
namespace mpl = boost::mpl;

template <int>
struct AbortConditionTraits;

template <>
struct AbortConditionTraits<surface>
{
  template <typename Result>
  static std::function<bool(const Result&)>
  create(const Surface& s)
  {
    return [s](const Result& r) -> bool {
      return r.endParameters.associatedSurface() == s;
    };
  }
};

template <typename Result>
class AbortCheck
{
  void
  addCondition(std::function<bool(Result&)> f)
  {
    m_conditions.push_back(std::move(f));
  }

private:
  std::vector<std::function<bool(Result&)>> m_conditions;
};

template <int Condition, typename R>
struct CreatorImpl;

template <int cond, typename R>
struct CreatorImpl<cond, R>
{
  static void
  apply(AbortCheck<R>& ac)
  {
    constexpr int lsb = cond & ~(cond - 1);
    ac.addCondition(AbortConditionTraits<lsb>::create());
    CreatorImpl<cond ^ lsb>::apply(ac);
  }
};

template <typename R>
struct CreatorImpl<0, R>
{
  static void
  apply(AbortCheck<R>&)
  {
  }
};

struct empty
{
};


template <template <typename...> class target,
          template <typename...> class sequence>
struct sequence_to_parameter_pack;

template <template <typename...> class target, typename... types>
struct sequence_to_parameter_pack<target, sequence<types...>>
{
  typedef typename target<types...> type;
};

template <int Condition, typename custom_result_extensions>
struct Creator
{
  typedef typename ResultTypeCreator<Condition, custom_result_extensions>::type
      result_type;
  static AbortCheck<R>
  create()
  {
    AbortCheck<R> result;
    CreatorImpl<Condition, R>::apply(result);
    return result;
  }
};
