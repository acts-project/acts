/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <functional>
#include <type_traits>
#include <utility>

namespace traccc {
template <typename T>
using rvalue_or_const_lvalue = std::disjunction<
    std::is_rvalue_reference<T>,
    std::conjunction<std::is_lvalue_reference<T>,
                     std::is_const<std::remove_reference_t<T>>>>;

template <typename T>
class algorithm {};

/**
 * @brief Unified algorithm semantics which convert an input to an output.
 *
 * This class provides a single, unified semantic for algorithms that can be
 * used throughout the traccc code. This virtual class is templated with an
 * input type and an output type, and any class implementing it is expected
 * to be able to transform the input type into the output type in one way or
 * another.
 *
 * @param A The input types of the algorithm
 * @param O The output type of the algorithm
 */
template <typename R, typename... A>
class algorithm<R(A...)> {
    public:
    using output_type = R;

    using function_type = R(A...);

    static_assert(
        std::conjunction<rvalue_or_const_lvalue<A>...>::value,
        "All arguments must be either affine types (rvalue references), or "
        "immutable constant types (const lvalue references).");

    virtual ~algorithm() = default;

    virtual output_type operator()(A... args) const = 0;
};

template <typename A, typename B, typename C, typename... R>
auto compose(
    std::function<std::remove_const_t<std::remove_reference_t<B>>(A)> f,
    std::function<C(B)> g, R... rs) {
    if constexpr (sizeof...(R) > 0) {
        auto h = compose(g, rs...);
        return [f, h](A&& i) { return h(f(std::forward<A>(i))); };
    } else {
        return [f, g](A&& i) { return g(f(std::forward<A>(i))); };
    }
}

template <typename A, typename B, typename C>
std::function<B(A&&)> side_effect(std::function<B(A)> f,
                                  std::function<void(const C&)> s) {
    return [=](A&& i) -> B {
        s(i);
        return f(std::forward<A>(i));
    };
}
}  // namespace traccc
