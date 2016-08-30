#ifndef ACTS_HAS_DUPLICATES_HPP
#define ACTS_HAS_DUPLICATES_HPP 1

namespace Acts {

namespace detail {

  template <typename... Args>
  struct has_duplicates;

  template <typename last>
  struct has_duplicates<last>
  {
    static constexpr bool value = false;
  };

  template <typename first, typename second, typename... others>
  struct has_duplicates<first, second, others...>
  {
  private:
    static constexpr bool _first  = has_duplicates<first, others...>::value;
    static constexpr bool _second = has_duplicates<second, others...>::value;

  public:
    static constexpr bool value = _first or _second;
  };

  template <typename first, typename... others>
  struct has_duplicates<first, first, others...>
  {
    static constexpr bool value = true;
  };

  template <typename... Args>
  constexpr bool has_duplicates_v = has_duplicates<Args...>::value;
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_HAS_DUPLICATES_HPP
