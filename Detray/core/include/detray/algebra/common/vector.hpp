// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

// @TODO: Remove this when Vc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#endif

// Project include(s)
#include "detray/algebra/common/array_operators.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace detray::algebra {

namespace storage {

namespace detail {
/// Make sure the vector/matrix dimension aligns with simd sizes
/// @see
/// https://gitlab.in2p3.fr/CodeursIntensifs/Fast5x5/-/blob/master/fast5x5.hpp?ref_type=heads
DETRAY_HOST_DEVICE
consteval std::size_t nearest_power_of_two(std::size_t min_value,
                                           std::size_t current_value) {
  // Computes the nearest power of two relative to `min_value` starting from
  // the power of two `current_value`
  return min_value <= current_value
             ? current_value
             : nearest_power_of_two(min_value, current_value * 2u);
}

}  // namespace detail

/// Vector wrapper for AoS vs interleaved SoA data. @c scalar_t can e.g. be a
/// SIMD vector.
template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
class DETRAY_ALIGN(
    alignof(array_t<scalar_t, detail::nearest_power_of_two(N, 2u)>)) vector {
 public:
  /// @returns the required size of the underlying array storage
  DETRAY_HOST_DEVICE
  static consteval std::size_t simd_size() {
    return concepts::value<scalar_t> ? detail::nearest_power_of_two(N, 2u) : N;
  }

  /// Scalar value type, or SIMD lane type for vectorized scalar storage.
  using scalar_type = scalar_t;
  /// Underlying data array type
  using array_type = array_t<scalar_t, simd_size()>;

  /// Default constructor sets all entries to zero.
  DETRAY_HOST_DEVICE
  constexpr vector() {
    if constexpr (!concepts::simd_scalar<scalar_type>) {
      zero_fill(std::make_index_sequence<simd_size()>{});
    }
  }

  /// Construct a vector from scalar elements.
  ///
  /// @{

  /// Construct vector in SoA layout from SIMD scalar elements.
  /// @param scals Element values, one for each logical vector component.
  template <typename... Scalars>
    requires(concepts::simd_scalar<scalar_t> && (sizeof...(Scalars) == N) &&
             ((concepts::simd_scalar<Scalars> ||
               std::convertible_to<Scalars, scalar_t>) &&
              ...))
  // Element-wise construction, not a conversion: the requires-clause pins the
  // argument count to N.
  // NOLINTNEXTLINE(google-explicit-constructor)
  DETRAY_HOST_DEVICE constexpr vector(Scalars &&...scals)
      : m_data{std::forward<Scalars>(scals)...} {}

  /// In order to avoid uninitialized values, which deteriorate the
  /// performance in explicitly vectorized code, the underlying data array is
  /// filled with zeroes if too few arguments are given.
  ///
  /// @param vals Element values; any unused storage lanes are initialized to
  /// zero where supported by the storage size.
  template <typename... Values>
    requires(!concepts::simd_scalar<scalar_t> && (sizeof...(Values) > 1) &&
             ((concepts::value<Values> ||
               std::convertible_to<Values, scalar_t>) &&
              ...))
  // NOLINTNEXTLINE(google-explicit-constructor)
  DETRAY_HOST_DEVICE constexpr vector(Values &&...vals) {
    static_assert(sizeof...(Values) <= N);

    // Fill up last entries, if necessary (explicitly for now)
    if constexpr ((simd_size() - N) == 1) {
      m_data = {std::forward<Values>(vals)..., 0.f};
    } else if constexpr ((simd_size() - N) == 2) {
      m_data = {std::forward<Values>(vals)..., 0.f, 0.f};
    } else if constexpr (sizeof...(Values) < simd_size()) {
      // @TODO: Does this actually work, yet?
      zero_fill(std::make_index_sequence<simd_size() - sizeof...(Values)>{});
    } else {
      m_data = {std::forward<Values>(vals)...};
    }
  }

  /// @}

  /// Construct from existing array storage.
  ///
  /// Implicit in both directions on purpose: together with the conversion
  /// operators below this makes the wrapper a drop-in for the array it wraps.
  /// @{
  /// @param vals Storage array to move into the wrapper.
  DETRAY_HOST_DEVICE
  // NOLINTNEXTLINE(google-explicit-constructor)
  constexpr vector(array_type &&vals) : m_data{std::move(vals)} {}

  /// @param vals Storage array to copy into the wrapper.
  DETRAY_HOST_DEVICE
  // NOLINTNEXTLINE(google-explicit-constructor)
  constexpr vector(const array_type &vals) : m_data{vals} {}
  /// @}

  /// Assignment operator from a vector with the same underlying storage.
  ///
  /// @param lhs Vector whose storage is assigned to this object.
  /// @returns A const reference to this vector after assignment.
  template <std::size_t M>
    requires(vector<N, scalar_t, array_t>::simd_size() ==
             vector<M, scalar_t, array_t>::simd_size())
  DETRAY_HOST_DEVICE constexpr const vector &operator=(
      const vector<M, scalar_t, array_t> &lhs) {
    m_data = lhs;
    return *this;
  }

  /// Convert the wrapper to its underlying mutable data array.
  /// @{
  /// @returns A mutable reference to the stored array.
  DETRAY_HOST_DEVICE
  // NOLINTNEXTLINE(google-explicit-constructor)
  constexpr operator array_type &() { return m_data; }
  /// Convert the wrapper to its underlying const data array.
  /// @returns A const reference to the stored array.
  DETRAY_HOST_DEVICE
  // NOLINTNEXTLINE(google-explicit-constructor)
  constexpr operator const array_type &() const { return m_data; }
  /// @}

  /// Return the underlying storage array.
  /// @returns A const reference to the stored array.
  DETRAY_HOST_DEVICE
  constexpr const auto &get() const { return m_data; }

  /// Access a stored component by index.
  /// @{
  /// @param i Zero-based component index.
  /// @returns A reference to the selected component.
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator[](std::size_t i) { return m_data[i]; }
  /// @param i Zero-based component index.
  /// @returns The selected component.
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator[](std::size_t i) const { return m_data[i]; }
  /// @}

  /// Multiply every component by a scalar factor.
  ///
  /// @param factor Scalar multiplier.
  /// @return The result of the underlying array multiplication.
  DETRAY_HOST_DEVICE
  constexpr decltype(auto) operator*=(scalar_type factor) noexcept {
    return m_data *= factor;
  }

  /// Equality operator
  /// @{
  /// AoS
  template <concepts::scalar S = scalar_t>
    requires(!concepts::simd_scalar<S>)
  DETRAY_HOST_DEVICE constexpr friend bool operator==(
      const vector &lhs, const vector &rhs) noexcept {
    const auto comp = lhs.compare(rhs);
    bool is_full = true;

    DETRAY_UNROLL_N(N)
    for (unsigned int i{0u}; i < N; ++i) {
      is_full = is_full && comp[i];
    }

    return is_full;
  }

  /// SoA
  template <concepts::scalar S = scalar_t>
    requires(concepts::simd_scalar<S>)
  DETRAY_HOST_DEVICE constexpr friend bool operator==(
      const vector &lhs, const vector &rhs) noexcept {
    const auto comp = lhs.compare(rhs);
    bool is_full = true;

    DETRAY_UNROLL_N(N)
    for (unsigned int i{0u}; i < N; ++i) {
      // Ducktyping the Vc::Vector::MaskType
      is_full = is_full && comp[i].isFull();
    }

    return is_full;
  }
  /// @}

  /// Compare this vector with another value for inequality.
  /// @param rhs Value to compare against.
  /// @returns Whether the values compare unequal.
  template <typename other_t>
  DETRAY_HOST_DEVICE constexpr bool operator!=(
      const other_t &rhs) const noexcept {
    return ((*this == rhs) == false);
  }

  /// Compare corresponding elements; SoA storage can produce a mask per lane.
  /// @param rhs Value whose elements are compared with this vector.
  /// @returns The per-element comparison results.
  template <typename other_t>
  DETRAY_HOST_DEVICE constexpr auto compare(const other_t &rhs) const noexcept {
    using result_t = decltype(m_data[0] == rhs[0]);

    std::array<result_t, N> comp;

    DETRAY_UNROLL_N(N)
    for (unsigned int i{0u}; i < N; ++i) {
      comp[i] = (m_data[i] == rhs[i]);
    }

    return comp;
  }

  /// Storage holding the value for every vector element.
  array_t<scalar_t, simd_size()> m_data;

 private:
  /// Sets the trailing uninitialized values to zero.
  template <std::size_t... Is>
  DETRAY_HOST_DEVICE constexpr void zero_fill(
      std::index_sequence<Is...> /*unused*/) noexcept {
    ((m_data[simd_size() - sizeof...(Is) + Is] = static_cast<scalar_t>(0)),
     ...);
  }
};

/// Macro declaring all instances of a specific arithmetic operator
#define DECLARE_VECTOR_OPERATORS(OP)                                           \
  template <std::size_t N, concepts::scalar scalar_t, concepts::value value_t, \
            template <typename, std::size_t> class array_t>                    \
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator OP(                     \
      const vector<N, scalar_t, array_t> &lhs, value_t rhs) noexcept {         \
    return lhs.m_data OP static_cast<scalar_t>(rhs);                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t, concepts::value value_t, \
            template <typename, std::size_t> class array_t>                    \
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator OP(                     \
      value_t lhs, const vector<N, scalar_t, array_t> &rhs) noexcept {         \
    return static_cast<scalar_t>(lhs) OP rhs.m_data;                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t>                    \
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator OP(                     \
      const vector<N, scalar_t, array_t> &lhs,                                 \
      const vector<N, scalar_t, array_t> &rhs) noexcept {                      \
    return lhs.m_data OP rhs.m_data;                                           \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t, typename other_t>  \
    requires(concepts::vector<other_t> || concepts::simd_scalar<other_t>)      \
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator OP(                     \
      const vector<N, scalar_t, array_t> &lhs, const other_t &rhs) noexcept {  \
    return lhs.m_data OP rhs;                                                  \
  }                                                                            \
  template <std::size_t N, concepts::scalar scalar_t,                          \
            template <typename, std::size_t> class array_t, typename other_t>  \
    requires(concepts::vector<other_t> || concepts::simd_scalar<other_t>)      \
  DETRAY_HOST_DEVICE constexpr decltype(auto) operator OP(                     \
      const other_t &lhs, const vector<N, scalar_t, array_t> &rhs) noexcept {  \
    return lhs OP rhs.m_data;                                                  \
  }

// Implement all arithmetic operations on top of @c vector.
// clang-format off
DECLARE_VECTOR_OPERATORS(+)
DECLARE_VECTOR_OPERATORS(-)
DECLARE_VECTOR_OPERATORS(*)
DECLARE_VECTOR_OPERATORS(/)
// clang-format on

// Clean up.
#undef DECLARE_VECTOR_OPERATORS

}  // namespace storage

namespace detail {

template <typename T>
struct is_storage_vector : public std::false_type {};

template <std::size_t N, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
struct is_storage_vector<algebra::storage::vector<N, scalar_t, array_t>>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_storage_vector_v = is_storage_vector<T>::value;

}  // namespace detail

}  // namespace detray::algebra
