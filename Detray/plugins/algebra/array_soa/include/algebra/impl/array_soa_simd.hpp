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
#include <cstddef>
#include <limits>
#include <ostream>
#include <type_traits>

namespace detray::algebra::array_soa {

/// Result of a lane-wise comparison: one boolean per lane
template <std::size_t W>
struct mask {
  static_assert(W > 0u, "A mask needs at least one lane");

  using value_type = bool;

  /// Holds the boolean for every lane
  std::array<bool, W> m_data{};

  /// Default constructor: all lanes false
  constexpr mask() = default;

  /// Broadcast constructor
  // NOLINTNEXTLINE(google-explicit-constructor)
  DETRAY_HOST_DEVICE constexpr mask(bool b) {
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      m_data[i] = b;
    }
  }

  /// Construct from the underlying array
  DETRAY_HOST_DEVICE constexpr explicit mask(const std::array<bool, W> &d)
      : m_data{d} {}

  /// @returns the number of lanes
  DETRAY_HOST_DEVICE
  static consteval std::size_t size() { return W; }

  /// Lane access
  /// @{
  DETRAY_HOST_DEVICE
  constexpr bool operator[](std::size_t i) const { return m_data[i]; }
  DETRAY_HOST_DEVICE
  constexpr bool &operator[](std::size_t i) { return m_data[i]; }
  /// @}

  /// @returns true if every lane is set
  DETRAY_HOST_DEVICE
  constexpr bool isFull() const {
    bool ret{true};
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret = ret && m_data[i];
    }
    return ret;
  }

  /// Logical operators (lane-wise, no short-circuit)
  /// @{
  DETRAY_HOST_DEVICE
  friend constexpr mask operator!(const mask &m) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = !m.m_data[i];
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator&&(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = a.m_data[i] && b.m_data[i];
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator||(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = a.m_data[i] || b.m_data[i];
    }
    return ret;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator&(const mask &a, const mask &b) {
    return a && b;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator|(const mask &a, const mask &b) {
    return a || b;
  }

  DETRAY_HOST_DEVICE
  friend constexpr mask operator^(const mask &a, const mask &b) {
    mask ret;
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret.m_data[i] = a.m_data[i] != b.m_data[i];
    }
    return ret;
  }
  /// @}

  /// Equality: true if all lanes agree
  /// @{
  DETRAY_HOST_DEVICE
  friend constexpr bool operator==(const mask &a, const mask &b) {
    bool ret{true};
    DETRAY_UNROLL_N(W)
    for (std::size_t i = 0u; i < W; ++i) {
      ret = ret && (a.m_data[i] == b.m_data[i]);
    }
    return ret;
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
      out << (m.m_data[i] ? "1" : "0");
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
template <concepts::value T, std::size_t W>
struct DETRAY_ALIGN(16) simd {
  static_assert(W > 0u, "A lane bundle needs at least one lane");

  using value_type = T;
  using mask_type = mask<W>;

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
      ret.m_data[i] = a.m_data[i] OP b.m_data[i];                        \
    }                                                                    \
    return ret;                                                          \
  }                                                                      \
  DETRAY_HOST_DEVICE                                                     \
  friend constexpr mask_type operator OP(const simd &a, T b) {           \
    mask_type ret;                                                       \
    DETRAY_UNROLL_N(W)                                                   \
    for (std::size_t i = 0u; i < W; ++i) {                               \
      ret.m_data[i] = a.m_data[i] OP b;                                  \
    }                                                                    \
    return ret;                                                          \
  }                                                                      \
  DETRAY_HOST_DEVICE                                                     \
  friend constexpr mask_type operator OP(T a, const simd &b) {           \
    mask_type ret;                                                       \
    DETRAY_UNROLL_N(W)                                                   \
    for (std::size_t i = 0u; i < W; ++i) {                               \
      ret.m_data[i] = a OP b.m_data[i];                                  \
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
