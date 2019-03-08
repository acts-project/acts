// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <string>
#include <type_traits>

namespace detail {

/// This file contains an implementation of the detection idiom in C++.
/// It's not currently in the standard, but can be implemented using
/// standard features. This implementation is largely taken from the C++
/// [technical specifications, library fundamentals v2(https://en.cppreference.com/w/cpp/experimental/is_detected)

/// Helper struct which cannot be constructe (or destroyes)d at all.
struct nonesuch
{
  ~nonesuch()               = delete;
  nonesuch(nonesuch const&) = delete;
  void
  operator=(nonesuch const&)
      = delete;
};

/**
 * The detector pattern works like this: there is a default type, that
 * accepts an "operation" that can be basically anything. It also accepts 
 * variadic arguments for that operation. The default type
 * has a member type that is std::false_type to indicate success or 
 * failure. It also has a member type "type" which captures a type result.
 * Then there is a specialization which attempts to instantiate the operation
 * with the given parameters, and tries to assign it into std::void_t. If the operation
 * fails to instantiate (say, the checked for type does not exist), the specialization
 * will not be instantiated, and the compiler falls back to the default type which contains
 * std::false_type. Since it happens inside while the compiler tries to find a better
 * matching template specialization than the default one (so basically overload resolution),
 * a compile error inside the operation is handled as a substitution failure, 
 * and is not an error. If the instantiation succeeds, the specialization contains a 
 * std::true_type, and an alias to the result of the operation.
 */

/**
 * This is the default specialization.
 * It does not attempt to instantiate `Op<Args...>` at all.
 *
 * 
 */
template <class Default,
          class AlwaysVoid,
          template <class...> class Op,
          class... Args>
struct detector
{
  using value_t = std::false_type;
  using type    = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...>
{
  // Note that std::void_t is a C++17 feature
  using value_t = std::true_type;
  using type    = Op<Args...>;
};

}  // namespace detail

template <template <class...> class Op, class... Args>
using is_detected =
    typename detail::detector<detail::nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t =
    typename detail::detector<detail::nonesuch, void, Op, Args...>::type;

template <class Expected, template <class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class Default, template <class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template <bool... Bs>
constexpr bool require = std::conjunction<std::bool_constant<Bs>...>::value;

template <template <class...> class Op, class... Args>
constexpr bool exists = is_detected<Op, Args...>::value;

template <class Exact, template <class...> class Op, class... Args>
constexpr bool identical_to = is_detected_exact<Exact, Op, Args...>::value;

template <typename T>
struct add_const
{
  using type = T const;
};

template <typename T>
using add_const_t = typename add_const<T>::type;

#define METHOD_TRAIT(trait_name, method_name)                                  \
  template <class T, typename R, typename... Arguments>                        \
  struct trait_name                                                            \
  {                                                                            \
    template <typename T_>                                                     \
    static constexpr bool is_const                                             \
        = not std::is_same_v<std::remove_const_t<T_>, T_>;                     \
                                                                               \
    template <typename T_, typename = int>                                     \
    struct fptr_meta                                                           \
    {                                                                          \
      template <typename... Arguments_>                                        \
      using type = typename std::                                              \
          integral_constant<decltype(std::declval<T_>().method_name(           \
                                std::declval<Arguments_>()...)) (T_::*)(       \
                                Arguments_...),                                \
                            &T_::method_name>::value_type;                     \
    };                                                                         \
                                                                               \
    template <typename T_>                                                     \
    struct fptr_meta<T_, std::enable_if_t<is_const<T_>, int>>                  \
    {                                                                          \
      template <typename... Arguments_>                                        \
      using type = typename std::                                              \
          integral_constant<decltype(std::declval<T_>().method_name(           \
                                std::declval<Arguments_>()...)) (T_::*)(       \
                                Arguments_...) const,                          \
                            &T_::method_name>::value_type;                     \
    };                                                                         \
                                                                               \
    template <typename T_, typename... Arguments_>                             \
    using fptr_meta_t = typename fptr_meta<T_>::template type<Arguments_...>;  \
    template <typename T_, typename... Arguments_>                             \
    using qual_ret = decltype(                                                 \
        std::declval<T_>().method_name(std::declval<Arguments_>()...));        \
                                                                               \
    template <typename T_, typename = int>                                     \
    struct tv                                                                  \
    {                                                                          \
      static constexpr bool value = false;                                     \
    };                                                                         \
    template <typename T_>                                                     \
    struct tv<T_,                                                              \
              std::enable_if_t<is_detected_exact<R,                            \
                                                 qual_ret,                     \
                                                 T_,                           \
                                                 Arguments...>::value,         \
                               int>>                                           \
    {                                                                          \
      static constexpr bool value                                              \
          = is_detected<fptr_meta_t, T, Arguments...>::value;                  \
    };                                                                         \
  }

template <typename T,
          typename R,
          template <class...> class M,
          typename... Arguments>
constexpr bool has_method = M<T, R, Arguments...>::template tv<T>::value;

#define MEMBER_TRAIT(trait_name, member_name)                                  \
  template <typename T>                                                        \
  using trait_name = decltype(std::declval<T>().member_name);

template <typename T, template <class...> class M, typename V>
constexpr bool has_member = identical_to<V, M, T>;

// -------
