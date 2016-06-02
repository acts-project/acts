#ifndef ACTS_GET_POSITION_H
#define ACTS_GET_POSITION_H 1

namespace Acts
{
  /// @cond detail
  namespace detail
  {
    /**
     * @brief get position of integral constant in template parameter pack
     *
     * @tparam T integral type of the values to be investigated
     * @tparam target target value whose position in the template parameter pack should be determined
     * @tparam values template parameter pack containing the list of values
     *
     * @return `get_position<T,target,values...>::value` yields the position of `target` inside `values`.
     *         If `target` is not in the list of `values`, a compile-time error is generated.
     */
    template<typename T, T target, T ... values>
    struct get_position;

    /// @cond
    template<typename T, T target, T ... others>
    struct get_position<T, target, target, others...>
    {
      enum
      {
        value = 0
      };
    };

    template<typename T, T target, T next, T ... others>
    struct get_position<T, target, next, others...>
    {
      enum
      {
        value = get_position<T, target, others...>::value + 1
      };
    };
    /// @endcond
  }  // end of namespace detail
  /// @endcond
}  // end of namespace Acts

#endif // ACTS_GET_POSITION_H
