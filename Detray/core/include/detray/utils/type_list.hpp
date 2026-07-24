// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/utils/string_helpers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <string>
#include <string_view>
#include <type_traits>

namespace detray {

namespace types {

/// @brief type list implementation
/// @see https://www.codingwiththomas.com/blog/getting-started-with-typelists
template <typename... Ts>
struct list {};

/// Number of types in the list
/// @{
template <typename = void>
struct get_size {};

template <typename... Ts>
struct get_size<list<Ts...>>
    : std::integral_constant<std::size_t, sizeof...(Ts)> {};

template <typename L>
constexpr inline std::size_t size{get_size<L>()};
/// @}

/// Does the list contain a particular type?
/// @{
template <typename T, std::size_t = 0, typename = void>
struct contains_impl : public std::false_type {};

template <typename U, std::size_t I, typename T, typename... Ts>
struct contains_impl<U, I, list<T, Ts...>>
    : public contains_impl<U, I + 1u, list<Ts...>> {};

template <typename T, std::size_t I, typename... Ts>
struct contains_impl<T, I, list<T, Ts...>> : public std::true_type {
  static constexpr std::size_t pos{I};
};

template <typename T, std::size_t I>
struct contains_impl<T, I, list<>> : public std::false_type {
  static constexpr auto pos{std::numeric_limits<std::size_t>::max()};
};

template <typename L, typename T>
inline constexpr bool contains = contains_impl<T, 0u, L>::value;

template <typename L, typename T>
inline constexpr std::size_t position = contains_impl<T, 0u, L>::pos;
/// @}

/// Access the first type
/// @{
template <typename = void>
struct get_front {};

template <typename T, typename... Ts>
struct get_front<list<T, Ts...>> {
  using type = T;
};
template <typename L>
using front = typename get_front<L>::type;
/// @}

/// Access the last type
/// @{
template <typename = void>
struct get_back {};

template <typename T, typename... Ts>
struct get_back<list<T, Ts...>> {
  using type = std::conditional_t<sizeof...(Ts) == 0, T,
                                  typename get_back<list<Ts...>>::type>;
};
// Base case
template <>
struct get_back<list<>> {
  using type = void;
};
template <typename L>
using back = typename get_back<L>::type;
/// @}

/// Access the N-th type
/// @{
template <std::size_t N, typename = void>
struct get_at {};

template <std::size_t N, typename T, typename... Ts>
struct get_at<N, list<T, Ts...>> {
  using type = std::remove_cvref_t<decltype(detray::get<N>(
      std::declval<dtuple<T, Ts...>>()))>;
};
template <typename L, std::size_t N>
using at = typename get_at<N, L>::type;
/// @}

/// Append a type
/// @{
template <typename N, typename = void>
struct do_push_back {};

template <typename N, typename... Ts>
struct do_push_back<N, list<Ts...>> {
  using type = list<Ts..., N>;
};
template <typename L, typename N>
using push_back = typename do_push_back<N, L>::type;

template <typename L, typename N>
using push_back_unique = std::conditional_t<contains<L, N>, L, push_back<L, N>>;
/// @}

/// Prepend a type
/// @{
template <typename N, typename = void>
struct do_push_front {};

template <typename N, typename... Ts>
struct do_push_front<N, list<Ts...>> {
  using type = list<N, Ts...>;
};
template <typename L, typename N>
using push_front = typename do_push_front<N, L>::type;
/// @}

/// Traits for the type list
/// @{
namespace detail {

template <typename = void>
struct is_type_list : public std::false_type {};

template <typename... Ts>
struct is_type_list<types::list<Ts...>> : public std::true_type {};

template <typename L>
inline constexpr bool is_type_list_v{is_type_list<L>::value};

}  // namespace detail
///@}

/// Print the type list
/// @{

/// @see
/// https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname
template <typename T>
std::string demangle_type_name() {
#if defined(__clang__)
  constexpr std::string_view prefix{"[T = "};
  constexpr std::string_view suffix{"]"};
  constexpr std::string_view function{__PRETTY_FUNCTION__};
#elif defined(__GNUC__)
  constexpr std::string_view prefix{"with T = "};
  constexpr std::string_view suffix{"; "};
  constexpr std::string_view function{__PRETTY_FUNCTION__};
#elif defined(_MSC_VER)
  constexpr std::string_view prefix{"get_type_name<"};
  constexpr std::string_view suffix{">(void)"};
  constexpr std::string_view function{__FUNCSIG__};
#else
  std::string type_name{typeid(T).name()};
  type_name += " ";
  std::string_view function{type_name};
  constexpr std::string_view prefix{""};
  constexpr std::string_view suffix{" "};
#endif

  const std::size_t start{function.find(prefix) + prefix.size()};
  const std::size_t end{function.find(suffix)};
  const std::size_t length{end - start};

  return std::string{function.substr(start, length)};
}

/// @returns the name of a type as string
/// @tparam T the type
template <typename T>
std::string get_name(bool full = false) {
  std::string tp_str{""};
  try {
    tp_str = detray::types::demangle_type_name<T>();
  } catch (...) {
    return "unknown";
  }

  if (tp_str.empty()) {
    return "unknown";
  }

  if (full) {
    return tp_str;
  }

  // Remove template parameter list by removing everything behind the first
  // occurrence of '<'
  auto tokens = detray::utils::split_at_delim(tp_str, '<');
  tp_str = tokens.front();

  tokens.clear();

  // Strip the namespaces by removing everything before the last occurrence of
  // ':'
  tokens = detray::utils::split_at_delim(tp_str, ':');
  return tokens.back();
}

template <typename = void>
struct print {
  print() { std::puts("Not a type_list!\n"); }
};

template <typename... Ts>
struct print<list<Ts...>> {
  template <typename P = void, typename... Ps>
  void print_typeid(bool full) {
    std::printf("%s", get_name<P>(full).c_str());

    // Keep unrolling the pack
    if constexpr (sizeof...(Ps) > 0) {
      std::puts(", ");
      return print_typeid<Ps...>(full);
    }
  }

  explicit print(bool full = true) {
    std::puts("type_list<");
    print_typeid<Ts...>(full);
    std::puts(">\n");
  }
};
/// @}

/// @brief Map the types of an existing type list to a new (smaller) type list
///
/// @tparam orig_list_t the original type list
/// @tparam type_selector type trait that maps a given type of the type list
///                       to a new type in a local typedef 'type'
/// @tparam I the current type position during recursion in the original list
/// @tparam Fs the types currently in the list (recursion)
///
/// @param [out] list a type list of selected types from the original list
///
/// @returns the filled type list
template <typename orig_list_t, class type_selector, std::size_t I = 0u,
          typename... Fs>
DETRAY_HOST_DEVICE consteval auto filtered_list(
    const list<Fs...>& filtered = list<>{}) {
  static_assert(sizeof...(Fs) <= I, "Can only map down to list with ");

  // The current list of mapped types
  using list_t = list<Fs...>;

  if constexpr (I == size<orig_list_t>) {
    return filtered;
  } else {
    // Map the current type corresponding to position 'I'
    // to the new type according to the given type_selector
    using next_type = at<orig_list_t, I>;
    using mapped_t = typename type_selector::template type<next_type>;

    // Don't add the same type multiple times
    if constexpr (contains<list_t, mapped_t>) {
      return filtered_list<orig_list_t, type_selector, I + 1u>(filtered);
    } else {
      return filtered_list<orig_list_t, type_selector, I + 1u>(
          push_back<list_t, mapped_t>{});
    }
  }
}

/// @brief Map the types of an existing type list to an array of indices
/// corresponding to the filtered type list (@see @c filtered_list )
///
/// @tparam orig_list_t the original type list
/// @tparam type_selector type trait that maps a given type of the type list
///                       to a new type in a local typedef 'type'
/// @tparam Fs the types currently in the list (recursion)
///
/// @param [out] list a type list of selected types from the original list
/// @param [out] idx_array a type list of selected types from the original list
///
/// @returns the filled index array
template <typename orig_list_t, class type_selector, std::size_t I = 0u,
          typename... Fs>
DETRAY_HOST_DEVICE consteval auto filtered_indices(
    const list<Fs...>& filtered,
    std::array<dindex, size<orig_list_t>> idx_array = {0}) {
  // The current list of mapped types
  using list_t = list<Fs...>;

  if constexpr (I == size<orig_list_t>) {
    return idx_array;
  } else {
    // Map the current type corresponding to position 'I'
    // to the new type according to the given type_selector
    using next_type = at<orig_list_t, I>;
    using mapped_t = typename type_selector::template type<next_type>;

    // If the type is already mapped, recover the index
    if constexpr (contains<list_t, mapped_t>) {
      idx_array[I] = position<list_t, mapped_t>;

      return filtered_indices<orig_list_t, type_selector, I + 1u>(filtered,
                                                                  idx_array);
    } else {
      idx_array[I] = sizeof...(Fs);

      return filtered_indices<orig_list_t, type_selector, I + 1u>(
          push_back<list_t, mapped_t>{}, idx_array);
    }
  }
}

}  // namespace types

/// Type list concepts
/// @{
namespace concepts {

template <typename L>
concept type_list = types::detail::is_type_list_v<L>;

}
///@}

}  // namespace detray
