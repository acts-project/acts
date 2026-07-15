/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::details::functor {
/**
 * @brief The identity functor, such that `identity<T>` is equivalent to `T`.
 */
template <typename T>
using identity = T;

/**
 * @brief The reference functor, such that `reference<T>` is equivalent to `T&`.
 */
template <typename T>
using reference = T&;

/**
 * @brief The const reference functor, such that `const_reference<T>` is
 * equivalent to `const T&`.
 */
template <typename T>
using const_reference = const T&;

/**
 * @brief A natural transformation between two functors.
 *
 * Given two functors F1 and F2 as well as some types Ts... such that
 * `F2<Ts...>` is a type, this produces `F1<Ts...>`.
 */
template <template <typename...> typename F, typename T>
struct reapply {};

template <template <typename...> typename F1,
          template <typename...> typename F2, typename... Ts>
struct reapply<F1, F2<Ts...>> {
    using type = F1<Ts...>;
};
}  // namespace traccc::details::functor
