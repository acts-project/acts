// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_concepts.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>

namespace detray::algebra::array_soa {

namespace detail {

/// Signed integer type with a size of @tparam S bytes
/// @{
template <std::size_t S>
struct lane_int {};

template <>
struct lane_int<1u> {
  using type = std::int8_t;
};

template <>
struct lane_int<2u> {
  using type = std::int16_t;
};

template <>
struct lane_int<4u> {
  using type = std::int32_t;
};

template <>
struct lane_int<8u> {
  using type = std::int64_t;
};
/// @}

/// Signed integer type with the same size as @tparam T
template <typename T>
using lane_int_t = typename lane_int<sizeof(T)>::type;

/// Unsigned integer type with the same size as @tparam T
template <typename T>
using lane_uint_t = std::make_unsigned_t<lane_int_t<T>>;

/// Largest alignment of a lane bundle in bytes: the size of the widest common
/// vector register (AVX-512) and of a cache line. A larger alignment does not
/// improve the vectorization, and the support of large extended alignments
/// depends on the compiler and the platform (e.g. on the stack in device code).
inline constexpr std::size_t max_lane_alignment{64u};

/// Alignment of @tparam W lanes of type @tparam T: The size of the lane bundle
/// rounded up to a power of two, so that one bundle fills a vector register.
/// The alignment is limited to @c max_lane_alignment. The size of a bundle is
/// still a multiple of the alignment, because both are powers of two.
template <typename T, std::size_t W>
inline constexpr std::size_t lane_alignment{std::bit_ceil(sizeof(T) * W) <
                                                    max_lane_alignment
                                                ? std::bit_ceil(sizeof(T) * W)
                                                : max_lane_alignment};

}  // namespace detail

/// Result of a lane-wise comparison: one integer per lane
///
/// A lane is set if all bits are set (-1) and not set if all bits are zero.
/// The lane integer has the same size as the values that were compared. This
/// is the same layout as a hardware vector mask, so that the compiler can
/// vectorize the comparisons, the logical operations and the masked
/// assignments.
///
/// @tparam W Number of lanes, must be a power of two.
/// @tparam I Signed integer type of a lane.
template <std::size_t W, std::signed_integral I>
struct DETRAY_ALIGN((detail::lane_alignment<I, W>)) mask {
  static_assert(std::has_single_bit(W),
                "The number of lanes of a mask must be a power of two");

  /// Value type of a single lane, when it is read
  using value_type = bool;
  /// Integer type that stores a lane
  using lane_type = I;

  /// Bit pattern of a lane that is set
  static constexpr I s_true{static_cast<I>(-1)};
  /// Bit pattern of a lane that is not set
  static constexpr I s_false{static_cast<I>(0)};

  /// Holds the bit pattern of every lane
  std::array<I, W> m_data{};

  /// Default constructor: all lanes false
  constexpr mask() = default;

  /// Broadcast constructor
  ///
  /// @param b Value for all lanes.
  // NOLINTNEXTLINE(google-explicit-constructor)
  DETRAY_HOST_DEVICE constexpr mask(bool b) {
    const I v{b ? s_true : s_false};
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] = v;
    }
  }

  /// Construct from one boolean per lane
  ///
  /// @param d Value for every lane.
  DETRAY_HOST_DEVICE constexpr explicit mask(const std::array<bool, W> &d) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] = d[i] ? s_true : s_false;
    }
  }

  /// Lane-wise conversion from a mask with a different lane integer size
  ///
  /// @tparam J Lane integer type of the other mask.
  ///
  /// @param o Mask to convert.
  template <std::signed_integral J>
    requires(!std::same_as<I, J>)
  DETRAY_HOST_DEVICE constexpr explicit mask(const mask<W, J> &o) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      // Conversion between signed integers keeps the values 0 and -1
      m_data[i] = static_cast<I>(o.m_data[i]);
    }
  }

  /// @returns the number of lanes
  DETRAY_HOST_DEVICE
  static consteval std::size_t size() { return W; }

  /// @param i Lane index.
  /// @returns true if lane @p i is set
  DETRAY_HOST_DEVICE
  constexpr bool operator[](std::size_t i) const {
    return m_data[i] != s_false;
  }

  /// Set a single lane
  ///
  /// @param i Lane index.
  /// @param b New value of the lane.
  DETRAY_HOST_DEVICE
  constexpr void set(std::size_t i, bool b) {
    m_data[i] = b ? s_true : s_false;
  }

  /// @returns true if every lane is set
  DETRAY_HOST_DEVICE
  constexpr bool isFull() const {
    I ret{s_true};
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret = static_cast<I>(ret & m_data[i]);
    }
    return ret != s_false;
  }

  /// Logical operators: lane-wise bit operations, no short-circuit
  /// @{
  DETRAY_HOST_DEVICE
  friend constexpr mask operator!(const mask &m) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = static_cast<I>(~m.m_data[i]);
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator&(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = static_cast<I>(a.m_data[i] & b.m_data[i]);
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator|(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = static_cast<I>(a.m_data[i] | b.m_data[i]);
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator^(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = static_cast<I>(a.m_data[i] ^ b.m_data[i]);
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator&&(const mask &a, const mask &b) {
    return a & b;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator||(const mask &a, const mask &b) {
    return a | b;
  }
  /// @}

  /// Equality: true if all lanes agree
  /// @{
  DETRAY_HOST_DEVICE
  friend constexpr bool operator==(const mask &a, const mask &b) {
    return (!(a ^ b)).isFull();
  }

  DETRAY_HOST_DEVICE
  friend constexpr bool operator!=(const mask &a, const mask &b) {
    return !(a == b);
  }
  /// @}

  /// Print the mask
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &out, const mask &m) {
    out << "[";
    for (std::size_t i = 0u; i < W; ++i) {
      out << (m[i] ? "1" : "0");
      if (i != W - 1u) {
        out << ", ";
      }
    }
    out << "]";
    return out;
  }
};

/// Lane bundle: one value per SIMD lane, stored in a @c std::array.
///
/// This is the scalar type of the SoA plugin. Every operation is a plain loop
/// over the lanes, so that the compiler can vectorize it.
///
/// @tparam T Value type of a lane.
/// @tparam W Number of lanes, must be a power of two.
template <concepts::value T, std::size_t W>
struct DETRAY_ALIGN((detail::lane_alignment<T, W>)) simd {
  static_assert(std::has_single_bit(W),
                "The number of lanes of a lane bundle must be a power of two");

  using value_type = T;
  using mask_type = mask<W, detail::lane_int_t<T>>;

  /// Holds the value for every lane
  std::array<T, W> m_data{};

  /// Default constructor: all lanes zero
  constexpr simd() = default;

  /// Broadcast constructor
  // NOLINTNEXTLINE(google-explicit-constructor)
  template <typename U>
    requires(std::convertible_to<U, T>)
  DETRAY_HOST_DEVICE constexpr simd(U v) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] = static_cast<T>(v);
    }
  }

  /// Construct from the underlying array
  DETRAY_HOST_DEVICE constexpr explicit simd(const std::array<T, W> &d)
      : m_data{d} {}

  /// Lane-wise conversion from a bundle of a different value type
  template <concepts::value U>
    requires(!std::same_as<U, T>)
  DETRAY_HOST_DEVICE constexpr explicit simd(const simd<U, W> &o) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] = static_cast<T>(o.m_data[i]);
    }
  }

  /// @returns the number of lanes
  DETRAY_HOST_DEVICE
  static consteval std::size_t size() { return W; }

  /// Lane access
  /// @{
  DETRAY_HOST_DEVICE
  constexpr T operator[](std::size_t i) const { return m_data[i]; }
  DETRAY_HOST_DEVICE
  constexpr T &operator[](std::size_t i) { return m_data[i]; }
  /// @}

  /// Compound assignment
  /// @{
  DETRAY_HOST_DEVICE
  constexpr simd &operator+=(const simd &rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] += rhs.m_data[i];
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator-=(const simd &rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] -= rhs.m_data[i];
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator*=(const simd &rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] *= rhs.m_data[i];
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator/=(const simd &rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] /= rhs.m_data[i];
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator+=(T rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] += rhs;
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator-=(T rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] -= rhs;
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator*=(T rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] *= rhs;
    }
    return *this;
  }

  DETRAY_HOST_DEVICE
  constexpr simd &operator/=(T rhs) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] /= rhs;
    }
    return *this;
  }
  /// @}

  /// Unary minus
  DETRAY_HOST_DEVICE
  friend constexpr simd operator-(const simd &a) {
    simd ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = -a.m_data[i];
    }
    return ret;
  }

/// Arithmetic operators: bundle-bundle, bundle-value and value-bundle
#define DETRAY_ARRAY_SOA_ARITHMETIC_OP(OP)                          \
  DETRAY_HOST_DEVICE                                                \
  friend constexpr simd operator OP(const simd &a, const simd &b) { \
    simd ret;                                                       \
    DETRAY_UNROLL_N(W)                                              \
    for (std::size_t i = 0u; i < W; ++i) {                          \
      ret.m_data[i] = a.m_data[i] OP b.m_data[i];                   \
    }                                                               \
    return ret;                                                     \
  }                                                                 \
  DETRAY_HOST_DEVICE                                                \
  friend constexpr simd operator OP(const simd &a, T b) {           \
    simd ret;                                                       \
    DETRAY_UNROLL_N(W)                                              \
    for (std::size_t i = 0u; i < W; ++i) {                          \
      ret.m_data[i] = a.m_data[i] OP b;                             \
    }                                                               \
    return ret;                                                     \
  }                                                                 \
  DETRAY_HOST_DEVICE                                                \
  friend constexpr simd operator OP(T a, const simd &b) {           \
    simd ret;                                                       \
    DETRAY_UNROLL_N(W)                                              \
    for (std::size_t i = 0u; i < W; ++i) {                          \
      ret.m_data[i] = a OP b.m_data[i];                             \
    }                                                               \
    return ret;                                                     \
  }

  // clang-format off
  DETRAY_ARRAY_SOA_ARITHMETIC_OP(+)
  DETRAY_ARRAY_SOA_ARITHMETIC_OP(-)
  DETRAY_ARRAY_SOA_ARITHMETIC_OP(*)
  DETRAY_ARRAY_SOA_ARITHMETIC_OP(/)
  // clang-format on

#undef DETRAY_ARRAY_SOA_ARITHMETIC_OP

/// Comparison operators: return a mask
#define DETRAY_ARRAY_SOA_COMPARISON_OP(OP)                               \
  DETRAY_HOST_DEVICE                                                     \
  friend constexpr mask_type operator OP(const simd &a, const simd &b) { \
    mask_type ret;                                                       \
    DETRAY_UNROLL_N(W)                                                   \
    for (std::size_t i = 0u; i < W; ++i) {                               \
      ret.m_data[i] = (a.m_data[i] OP b.m_data[i]) ? mask_type::s_true   \
                                                   : mask_type::s_false; \
    }                                                                    \
    return ret;                                                          \
  }                                                                      \
  DETRAY_HOST_DEVICE                                                     \
  friend constexpr mask_type operator OP(const simd &a, T b) {           \
    mask_type ret;                                                       \
    DETRAY_UNROLL_N(W)                                                   \
    for (std::size_t i = 0u; i < W; ++i) {                               \
      ret.m_data[i] =                                                    \
          (a.m_data[i] OP b) ? mask_type::s_true : mask_type::s_false;   \
    }                                                                    \
    return ret;                                                          \
  }                                                                      \
  DETRAY_HOST_DEVICE                                                     \
  friend constexpr mask_type operator OP(T a, const simd &b) {           \
    mask_type ret;                                                       \
    DETRAY_UNROLL_N(W)                                                   \
    for (std::size_t i = 0u; i < W; ++i) {                               \
      ret.m_data[i] =                                                    \
          (a OP b.m_data[i]) ? mask_type::s_true : mask_type::s_false;   \
    }                                                                    \
    return ret;                                                          \
  }

  // clang-format off
  DETRAY_ARRAY_SOA_COMPARISON_OP(==)
  DETRAY_ARRAY_SOA_COMPARISON_OP(!=)
  DETRAY_ARRAY_SOA_COMPARISON_OP(<)
  DETRAY_ARRAY_SOA_COMPARISON_OP(<=)
  DETRAY_ARRAY_SOA_COMPARISON_OP(>)
  DETRAY_ARRAY_SOA_COMPARISON_OP(>=)
  // clang-format on

#undef DETRAY_ARRAY_SOA_COMPARISON_OP

  /// Print the bundle
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &out, const simd &v) {
    out << "[";
    for (std::size_t i = 0u; i < W; ++i) {
      out << v.m_data[i];
      if (i != W - 1u) {
        out << ", ";
      }
    }
    out << "]";
    return out;
  }
};

}  // namespace detray::algebra::array_soa

namespace std {

/// Numeric limits of a lane bundle: the constants of the value type, the
/// values broadcast to every lane
template <detray::concepts::value T, std::size_t W>
struct numeric_limits<detray::algebra::array_soa::simd<T, W>>
    : public numeric_limits<T> {
  using simd_t = detray::algebra::array_soa::simd<T, W>;

  DETRAY_HOST_DEVICE
  static constexpr simd_t min() noexcept {
    return simd_t{numeric_limits<T>::min()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t max() noexcept {
    return simd_t{numeric_limits<T>::max()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t lowest() noexcept {
    return simd_t{numeric_limits<T>::lowest()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t epsilon() noexcept {
    return simd_t{numeric_limits<T>::epsilon()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t round_error() noexcept {
    return simd_t{numeric_limits<T>::round_error()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t infinity() noexcept {
    return simd_t{numeric_limits<T>::infinity()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t quiet_NaN() noexcept {
    return simd_t{numeric_limits<T>::quiet_NaN()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t signaling_NaN() noexcept {
    return simd_t{numeric_limits<T>::signaling_NaN()};
  }
  DETRAY_HOST_DEVICE
  static constexpr simd_t denorm_min() noexcept {
    return simd_t{numeric_limits<T>::denorm_min()};
  }
};

}  // namespace std
