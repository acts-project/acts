// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/type_list.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

namespace detray::detail {

/// get function accessor for std::tuple
///
/// usage example:
/// detail::get<0>(tuple)
/// @{
using std::get;

template <std::size_t I, typename... value_types>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const ::detray::dtuple<value_types...>& tuple) noexcept {
  return ::detray::get<I>(tuple);
}

template <std::size_t I, typename... value_types>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    ::detray::dtuple<value_types...>& tuple) noexcept {
  return ::detray::get<I>(tuple);
}

template <typename query_t, typename... value_types>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const ::detray::dtuple<value_types...>& tuple) noexcept {
  return ::detray::get<get_type_pos_v<query_t, value_types...>>(tuple);
}

template <typename query_t, typename... value_types>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    ::detray::dtuple<value_types...>& tuple) noexcept {
  return ::detray::get<get_type_pos_v<query_t, value_types...>>(tuple);
}
/// @}

/// tuple_element for std::tuple
///
/// usage example:
/// detail::tuple_element< int, tuple_t >::type
/// @{
template <std::size_t N, class T>
struct tuple_element;

// std::tuple
template <std::size_t N, typename... value_types>
struct tuple_element<N, std::tuple<value_types...>>
    : std::tuple_element<N, std::tuple<value_types...>> {};

// detray::dtuple
/// Select tuple element `N` from a detray tuple.
template <std::size_t N, typename... value_types>
struct tuple_element<N, detray::dtuple<value_types...>> {
  /// Decayed type of the selected element.
  using type = std::decay_t<decltype(::detray::get<N>(
      std::declval<detray::dtuple<value_types...>>()))>;
};

template <std::size_t N, class T>
using tuple_element_t = typename tuple_element<N, T>::type;
/// @}

/// tuple_size for std::tuple
///
/// usage example:
/// detail::tuple_size< tuple_t >::value
/// @{
template <class T>
struct tuple_size;

// std::tuple
template <typename... value_types>
struct tuple_size<std::tuple<value_types...>>
    : std::tuple_size<std::tuple<value_types...>> {};

// detray::dtuple
/// Report the number of elements in a detray tuple.
template <typename... value_types>
struct tuple_size<::detray::dtuple<value_types...>> {
  /// Number of element types in the tuple.
  static constexpr std::size_t value = sizeof...(value_types);
};

template <class T>
constexpr std::size_t tuple_size_v{tuple_size<T>::value};
/// @}

/// make_tuple for std::tuple
/// users have to specify tuple_t for detail::make_tuple
///
/// usage example
/// detail::make_tuple<tuple_t>(args...)
/// @{
template <class T>
struct unwrap_refwrapper {
  /// Decayed type when the argument is not a reference wrapper.
  using type = T;
};

template <class T>
struct unwrap_refwrapper<std::reference_wrapper<T>> {
  /// Referenced type, preserving the wrapper's reference semantics.
  using type = T&;
};

template <class T>
using unwrap_decay_t = typename unwrap_refwrapper<std::decay_t<T>>::type;

// make_tuple for std::tuple
template <template <typename...> class tuple_t, class... value_types>
  requires std::is_same_v<tuple_t<value_types...>, std::tuple<value_types...>>
DETRAY_HOST constexpr std::tuple<unwrap_decay_t<value_types>...> make_tuple(
    value_types&&... args) {
  return std::make_tuple(std::forward<value_types>(args)...);
}

// make_tuple for detray::dtuple
template <template <typename...> class tuple_t, class... value_types>
  requires std::is_same_v<tuple_t<value_types...>,
                          detray::dtuple<value_types...>>
DETRAY_HOST_DEVICE constexpr detray::dtuple<unwrap_decay_t<value_types>...>
make_tuple(value_types&&... args) {
  return detray::dtuple<unwrap_decay_t<value_types>...>{
      std::forward<value_types>(args)...};
}
/// @}

/// Check if the tuple contains a type
/// @see
/// https://stackoverflow.com/questions/25958259/how-do-i-find-out-if-a-tuple-contains-a-type
/// @{
template <typename T, typename tuple_t>
struct has_type;

// std::tuple
template <typename T>
struct has_type<T, std::tuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

template <typename T, typename... Ts>
struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

// detray::dtuple
template <typename T>
struct has_type<T, detray::dtuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct has_type<T, detray::dtuple<U, Ts...>>
    : has_type<T, detray::dtuple<Ts...>> {};

template <typename T, typename... Ts>
struct has_type<T, detray::dtuple<T, Ts...>> : std::true_type {};

template <typename T, class tuple_t>
constexpr bool has_type_v = has_type<T, tuple_t>::value;
///@}

/// Concatenate tuple types
/// @{
template <typename... tuple_ts>
struct tuple_cat_type {};

template <typename... Args>
struct tuple_cat_type<std::tuple<Args...>> {
  using type = std::tuple<Args...>;
};

template <typename... Args1, typename... Args2, typename... tuple_ts>
struct tuple_cat_type<std::tuple<Args1...>, std::tuple<Args2...>, tuple_ts...> {
  using type = typename tuple_cat_type<std::tuple<Args1..., Args2...>,
                                       tuple_ts...>::type;
};

/// Preserve a detray tuple type as the concatenation result.
template <typename... Args>
struct tuple_cat_type<detray::dtuple<Args...>> {
  /// Tuple type produced by concatenation.
  using type = detray::dtuple<Args...>;
};

/// Concatenate two or more detray tuple type lists recursively.
template <typename... Args1, typename... Args2, typename... tuple_ts>
struct tuple_cat_type<detray::dtuple<Args1...>, detray::dtuple<Args2...>,
                      tuple_ts...> {
  /// Tuple type containing the concatenated element types.
  using type = typename tuple_cat_type<detray::dtuple<Args1..., Args2...>,
                                       tuple_ts...>::type;
};

template <typename... tuple_ts>
using tuple_cat_t = typename tuple_cat_type<tuple_ts...>::type;
/// @}

/// Remove duplicate types from tuple
/// @{
/// Convert tuple element types into a type list.
template <typename tuple_t>
struct tuple_to_type_list {};

template <template <typename...> class tuple_t, typename... Args>
struct tuple_to_type_list<tuple_t<Args...>> {
  /// Type list preserving the input tuple order.
  using type = detray::types::list<Args...>;
};

/// Convert a type list into the supplied tuple template.
template <template <typename...> class tuple_t, typename type_list_t>
struct type_list_to_tuple {};

template <template <typename...> class tuple_t, typename... Args>
struct type_list_to_tuple<tuple_t, detray::types::list<Args...>> {
  /// Tuple type preserving the input list order.
  using type = tuple_t<Args...>;
};

/// Collect the first occurrence of each input type in order.
template <typename result_list_t, typename input_list_t>
struct unique_types {};

/// Return the accumulated list after all input types have been visited.
template <typename result_list_t>
struct unique_types<result_list_t, detray::types::list<>> {
  /// Type list with repeated types removed.
  using type = result_list_t;
};

/// Append a previously unseen type and continue with the remaining inputs.
template <typename result_list_t, typename arg_t, typename... input_args>
struct unique_types<result_list_t, detray::types::list<arg_t, input_args...>> {
 private:
  using next_result_list_t =
      detray::types::push_back_unique<result_list_t, arg_t>;

 public:
  /// Type list with stable left-to-right uniqueness.
  using type = typename unique_types<next_result_list_t,
                                     detray::types::list<input_args...>>::type;
};

/// Remove repeated tuple types while retaining their first occurrence.
template <typename tuple_t>
struct unique_type {};

/// Preserve the original tuple template when collecting unique types.
template <template <typename...> class tuple_t, typename... Args>
struct unique_type<tuple_t<Args...>> {
 private:
  using unique_list_t = typename unique_types<
      detray::types::list<>,
      typename tuple_to_type_list<tuple_t<Args...>>::type>::type;

 public:
  /// Tuple type with stable left-to-right uniqueness.
  using type = typename type_list_to_tuple<tuple_t, unique_list_t>::type;
};

template <typename tuple_t>
using unique_t = typename unique_type<tuple_t>::type;
/// @}

/// Check for equality of tuple types modulo permutation
/// @{
template <typename T1, typename T2>
struct is_permutation : public std::false_type {};

template <>
struct is_permutation<detray::dtuple<>, detray::dtuple<>>
    : public std::true_type {};

/// Compare the two detray tuple type lists as unordered sets of types.
template <typename... Args1, typename... Args2>
struct is_permutation<detray::dtuple<Args1...>, detray::dtuple<Args2...>> {
  /// First tuple type being compared.
  using T1 = detray::dtuple<Args1...>;
  /// Second tuple type being compared.
  using T2 = detray::dtuple<Args2...>;

  /// Check whether every type in the pack occurs in the other tuple.
  /// @returns Whether the tuple types form the same set.
  template <typename o_tuple_t, typename T, typename... U>
  static consteval bool compare() {
    if constexpr (detray::detail::has_type_v<T, o_tuple_t>) {
      if constexpr (sizeof...(U) == 0u) {
        return true;
      } else {
        return compare<o_tuple_t, U...>();
      }
    } else {
      return false;
    }
  }

  /// Whether both tuples contain the same types, regardless of order.
  static constexpr bool value =
      ((detray::detail::tuple_size_v<T1> == detray::detail::tuple_size_v<T2>) &&
       (compare<T1, Args2...>() && compare<T2, Args1...>()));
};

template <typename T1, typename T2>
constexpr bool is_permutation_v = is_permutation<T1, T2>::value;
/// @}

/// Check trait on ever element of the tuple
/// @{

// Any
template <template <typename...> class trait, typename tuple_t>
struct tuple_any {};

template <template <typename...> class trait, typename... Args>
struct tuple_any<trait, std::tuple<Args...>> {
  static constexpr bool value = (trait<Args>::value || ...);
};
/// Apply `trait` to each element type of a detray tuple and combine with OR.
template <template <typename...> class trait, typename... Args>
struct tuple_any<trait, detray::dtuple<Args...>> {
  /// Whether at least one element type satisfies `trait`.
  static constexpr bool value = (trait<Args>::value || ...);
};

template <template <typename...> class trait, typename tuple_t>
constexpr bool tuple_any_v = tuple_any<trait, tuple_t>::value;

// All
template <template <typename...> class trait, typename tuple_t>
struct tuple_all {};

template <template <typename...> class trait, typename... Args>
struct tuple_all<trait, std::tuple<Args...>> {
  static constexpr bool value = (trait<Args>::value && ...);
};
/// Apply `trait` to each element type of a detray tuple and combine with AND.
template <template <typename...> class trait, typename... Args>
struct tuple_all<trait, detray::dtuple<Args...>> {
  /// Whether every element type satisfies `trait`.
  static constexpr bool value = (trait<Args>::value && ...);
};

template <template <typename...> class trait, typename tuple_t>
constexpr bool tuple_all_v = tuple_all<trait, tuple_t>::value;

/// @}

}  // namespace detray::detail
