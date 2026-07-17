/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>
#include <utility>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

template <typename... Ts>
struct tuple {};

template <typename T, typename... Ts>
struct tuple<T, Ts...> {

    TRACCC_HOST_DEVICE constexpr tuple(){};

    constexpr tuple(const tuple &o)
        requires(std::is_copy_constructible_v<T> &&
                 (std::is_copy_constructible_v<Ts> && ...))
    = default;

    template <typename U, typename... Us>
        requires std::is_constructible_v<T, U &&> &&
                     std::is_constructible_v<tuple<Ts...>, Us &&...>
    TRACCC_HOST_DEVICE explicit constexpr tuple(const tuple<U, Us...> &o)
        : v(o.v), r(o.r) {}

    constexpr tuple(tuple &&o) noexcept
        requires(std::is_move_constructible_v<T> &&
                 (std::is_move_constructible_v<Ts> && ...))
    = default;

    template <typename U, typename... Us>
        requires std::is_constructible_v<T, U &&> &&
                     std::is_constructible_v<tuple<Ts...>, Us &&...>
    TRACCC_HOST_DEVICE explicit constexpr tuple(tuple<U, Us...> &&o)
        : v(std::move(o.v)), r(std::move(o.r)) {}

    template <typename U, typename... Us>
        requires std::is_constructible_v<T, U &&> &&
                     std::is_constructible_v<tuple<Ts...>, Us &&...> &&
                     (!(std::is_same_v<tuple, U> ||
                        (std::is_same_v<tuple, Us> || ...)))
    TRACCC_HOST_DEVICE explicit constexpr tuple(U &&_v, Us &&..._r)
        : v(std::forward<U>(_v)), r(std::forward<Us>(_r)...) {}

    constexpr ~tuple() noexcept = default;

    constexpr tuple &operator=(const tuple &other)
        requires(std::is_copy_assignable_v<T> &&
                 (std::is_copy_assignable_v<Ts> && ...))
    = default;

    template <typename U, typename... Us>
    TRACCC_HOST_DEVICE constexpr tuple &operator=(const tuple<U, Us...> &other)
        requires(std::is_assignable_v<T &, const U &> &&
                 (std::is_assignable_v<Ts &, const Us &> && ...))
    {
        v = other.v;
        r = other.r;
        return *this;
    }

    constexpr tuple &operator=(tuple &&other) noexcept
        requires(std::is_move_assignable_v<T> &&
                 (std::is_move_assignable_v<Ts> && ...))
    = default;

    template <typename U, typename... Us>
    TRACCC_HOST_DEVICE constexpr tuple &operator=(tuple<U, Us...> &&other)
        requires(std::is_assignable_v<T &, U> &&
                 (std::is_assignable_v<Ts &, Us> && ...))
    {
        v = std::move(other.v);
        r = std::move(other.r);
        return *this;
    }

    T v;
    tuple<Ts...> r;
};
}  // namespace traccc
