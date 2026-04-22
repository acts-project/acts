// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/bit_encoder.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <array>
#include <concepts>
#include <cstdint>
#include <ostream>
#include <type_traits>

namespace detray {

template <typename I>
DETRAY_HOST inline std::ostream& operator<<(std::ostream& os,
                                            const darray<I, 2>& r) {
  os << "[0: " << r[0] << ", 1: " << r[1] << "]";
  return os;
}

namespace detail {

/// Weather to interpret index ranges as containing the lower index and the
/// number of elements or the lower and upper index of the range
inline constexpr bool sized_index_range{true};

/// @brief Index range
///
/// @tparam index_t the type of indices
/// @tparam value_t underlying bit encoded value type
/// @tparam lower_mask how many bist to use for the lower index
template <typename index_t, bool contains_size = !sized_index_range,
          typename value_t = std::uint_least64_t,
          value_t lower_mask = 0xffffffff00000000,
          value_t upper_mask = ~lower_mask>
  requires std::convertible_to<index_t, value_t>
struct index_range {
  using index_type = index_t;
  using encoder = detail::bit_encoder<value_t>;

  constexpr index_range() { set_zeros(); }

  /// Construct around encoded value @param v
  DETRAY_HOST_DEVICE
  explicit constexpr index_range(value_t v) : m_value{v} { set_zeros(); }

  /// Construct from a lower and an upper index/number of elements
  DETRAY_HOST_DEVICE
  constexpr index_range(index_t l, index_t u) {
    encoder::template set_bits<lower_mask>(m_value, static_cast<value_t>(l));
    encoder::template set_bits<upper_mask>(m_value, static_cast<value_t>(u));
    set_zeros();
  }

  /// Explicit conversion to underlying value type
  explicit constexpr operator value_t() const { return m_value; }

  /// Elementwise access.
  /// @{
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](const std::size_t i) {
    assert(i <= 1u);
    return (i == 0u) ? lower() : upper();
  }
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](const std::size_t i) const {
    assert(i <= 1u);
    return (i == 0u) ? lower() : upper();
  }
  /// @}

  /// @returns the size of the index range
  DETRAY_HOST_DEVICE
  constexpr auto size() const -> std::size_t {
    if constexpr (contains_size) {
      return static_cast<std::size_t>(
          encoder::template get_bits<upper_mask>(m_value));
    } else {
      assert(upper() >= lower());
      return static_cast<std::size_t>(upper() - lower());
    }
  }

  /// @return the first index of the range - const
  DETRAY_HOST_DEVICE
  constexpr auto lower() const -> index_type {
    return static_cast<index_type>(
        encoder::template get_bits<lower_mask>(m_value));
  }

  /// @return the last index of the range - const
  DETRAY_HOST_DEVICE
  constexpr auto upper() const -> index_type {
    const auto upper_val{static_cast<index_type>(
        encoder::template get_bits<upper_mask>(m_value))};

    if constexpr (contains_size) {
      return lower() + upper_val;
    } else {
      return upper_val;
    }
  }

  /// Set the lower index.
  DETRAY_HOST_DEVICE
  constexpr index_range& set_lower(index_type l) {
    encoder::template set_bits<lower_mask>(m_value, static_cast<value_t>(l));
    return *this;
  }

  /// Set the upper index.
  DETRAY_HOST_DEVICE
  constexpr index_range& set_upper(index_type u) {
    if constexpr (contains_size) {
      assert((u >= lower()) &&
             "Sized index range: Upper index can only be set if lower "
             "index is already set");

      // Set invalid size correctly
      constexpr auto inv_v{~static_cast<value_t>(0)};
      constexpr auto inv_upper{static_cast<index_type>(
          encoder::template get_bits<upper_mask>(inv_v))};
      const bool inv_idx{u == inv_upper};

      encoder::template set_bits<upper_mask>(
          m_value, static_cast<value_t>(inv_idx ? u : u - lower()));
    } else {
      encoder::template set_bits<upper_mask>(m_value, static_cast<value_t>(u));
    }
    return *this;
  }

  /// Shift both indices by @param s. Allow negative shift
  DETRAY_HOST_DEVICE
  constexpr index_range& shift(long s) {
    assert(static_cast<long>(lower()) + s >= 0l);

    // Update lower index
    encoder::template set_bits<lower_mask>(
        m_value, static_cast<value_t>(static_cast<long>(lower()) + s));
    // Shift upper index, if range contains upper index instead of size
    if constexpr (!contains_size) {
      encoder::template set_bits<upper_mask>(
          m_value, static_cast<value_t>(static_cast<long>(upper()) + s));
    }

    return *this;
  }

  /// Check whether the index range is valid to use.
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid() const noexcept {
    return encoder::template is_invalid<lower_mask, upper_mask>(m_value);
  }

  /// Comparison operators
  /// @{
  /// Equality operator
  constexpr bool operator==(const index_range& other) const {
    return (lower() == other.lower()) && (upper() == other.upper());
  };

  /// Comparison operators
  DETRAY_HOST_DEVICE
  friend constexpr auto operator<=>(const index_range lhs,
                                    const index_range rhs) noexcept {
    const auto l{lhs.size()};
    const auto r{rhs.size()};
    if (l < r || (l == r && l < r)) {
      return std::strong_ordering::less;
    }
    if (l > r || (l == r && l > r)) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equivalent;
  }
  /// @}

  /// Arithmetic operators
  /// @{
  /// @returns extended index range that contains both initial ranges
  DETRAY_HOST_DEVICE
  friend index_range operator+(const index_range lhs, const index_range rhs) {
    const index_t new_lower{lhs.lower() < rhs.lower() ? lhs.lower()
                                                      : rhs.lower()};
    const index_t new_upper{lhs.upper() > rhs.upper() ? lhs.upper()
                                                      : rhs.upper()};

    return new_lower < new_upper
               ? index_range{}.set_lower(new_lower).set_upper(new_upper)
               : index_range{}.set_lower(new_upper).set_upper(new_lower);
  }

  /// Expands upper range boundary by @param index
  DETRAY_HOST_DEVICE
  friend index_range operator+(const index_range lhs, const index_type index) {
    return index_range{}.set_lower(lhs.lower()).set_upper(lhs.upper() + index);
  }

  /// @returns index range that was reduced to minimal coverage
  DETRAY_HOST_DEVICE
  friend index_range operator-(const index_range lhs, const index_range rhs) {
    const index_t new_lower{lhs.lower() > rhs.lower() ? lhs.lower()
                                                      : rhs.lower()};
    const index_t new_upper{lhs.upper() < rhs.upper() ? lhs.upper()
                                                      : rhs.upper()};

    return new_lower < new_upper
               ? index_range{}.set_lower(new_lower).set_upper(new_upper)
               : index_range{}.set_lower(new_upper).set_upper(new_lower);
  }

  /// Shrinks upper range boundary by @param index
  DETRAY_HOST_DEVICE
  friend index_range operator-(const index_range lhs, const index_type& index) {
    const index_t new_lower{lhs.lower()};
    const index_t new_upper{lhs.upper() - index};

    return new_lower < new_upper
               ? index_range{}.set_lower(new_lower).set_upper(new_upper)
               : index_range{}.set_lower(new_upper).set_upper(new_lower);
  }
  /// @}

 private:
  /// Set bits that are not covered by the masks to zero, so that the
  /// assertion in the bit_encoder works (does nothing if the lower and upper
  /// masks cover all of the bits in the value_t)
  constexpr void set_zeros() {
    constexpr value_t uncovered_mask{~(lower_mask | upper_mask)};
    if constexpr (uncovered_mask != 0u) {
      encoder::template set_bits<uncovered_mask>(m_value, 0u);
    }
  }

  /// Print operator
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const index_range& ir) {
    if (ir.is_invalid()) {
      return (os << "undefined");
    }

    os << "[" << ir.lower() << ", " << ir.upper() << "]";
    return os;
  }

  /// The encoded value. Default: All bits set to 1 (invalid)
  value_t m_value = ~static_cast<value_t>(0);
};

/// @brief Simple multi-index structure
///
/// @tparam N number of indices that are held by this type
/// @tparam index_t type of indices
template <typename index_t, std::size_t N>
struct multi_index {
  using index_type = index_t;

  DETRAY_HOST_DEVICE
  constexpr multi_index() : m_indices{get_invalid_array()} {}

  /// Construct from a list of values
  template <typename... I>
    requires(sizeof...(I) == N &&
             std::conjunction_v<std::is_convertible<I, index_t>...>)
  DETRAY_HOST_DEVICE constexpr multi_index(I... indices)
      : m_indices{indices...} {}

  /// @returns the number of contained indices
  DETRAY_HOST_DEVICE
  constexpr static auto size() -> std::size_t { return N; }

  /// Check whether the multi index is valid to use.
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid() const noexcept {
    bool is_inv_elem{false};
    for (std::size_t i = 0u; i < size(); ++i) {
      is_inv_elem |= detail::is_invalid_value(m_indices[i]);
    }
    return is_inv_elem;
  }

  /// Elementwise access.
  /// @{
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](const std::size_t i) { return m_indices[i]; }
  DETRAY_HOST_DEVICE
  decltype(auto) operator[](const std::size_t i) const { return m_indices[i]; }
  /// @}

  /// Equality operator @returns true if all bin indices match.
  bool operator==(const multi_index& rhs) const = default;

 private:
  /// Compile time generate an array of invalid values
  static consteval darray<index_t, N> get_invalid_array() {
    darray<index_t, N> result{};
    for (std::size_t i = 0u; i < N; ++i) {
      result[i] = detail::invalid_value<index_t>();
    }
    return result;
  }

  /// Print operator
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const multi_index& mi) {
    os << "[";
    bool writeSeparator = false;
    for (auto i = 0u; i < N; ++i) {
      if (writeSeparator) {
        os << ", ";
      }
      os << i << ": " << mi[i];
      writeSeparator = true;
    }
    os << "]";

    return os;
  }

  /// The underlying index array
  darray<index_t, N> m_indices;
};

/// @brief Ties an object type and an index into a container together.
///
/// @tparam id_type Represents the type of object that is being indexed
/// @tparam index_type The type of indexing needed for the indexed type's
///         container (e.g. single index, range, multiindex)
template <typename id_t, typename index_t,
          typename value_t = std::uint_least32_t, value_t id_mask = 0xf0000000,
          value_t index_mask = ~id_mask>
  requires(std::convertible_to<value_t, index_t> ||
           std::constructible_from<index_t, value_t>)
struct typed_index {
  using id_type = id_t;
  using index_type = index_t;
  using encoder = detail::bit_encoder<value_t>;

  constexpr typed_index() = default;

  DETRAY_HOST_DEVICE
  constexpr typed_index(const id_t id, const index_t idx) {
    encoder::template set_bits<id_mask>(m_value, static_cast<value_t>(id));
    encoder::template set_bits<index_mask>(m_value, static_cast<value_t>(idx));
  }

  /// @return the type id - const
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> id_type {
    return static_cast<id_type>(encoder::template get_bits<id_mask>(m_value));
  }

  /// @return the index - const
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> index_type {
    return static_cast<index_type>(
        encoder::template get_bits<index_mask>(m_value));
  }

  /// Set the link id.
  DETRAY_HOST_DEVICE
  constexpr typed_index& set_id(id_type id) {
    encoder::template set_bits<id_mask>(m_value, static_cast<value_t>(id));
    return *this;
  }

  /// Set the link index.
  DETRAY_HOST_DEVICE
  constexpr typed_index& set_index(index_type index) {
    encoder::template set_bits<index_mask>(m_value,
                                           static_cast<value_t>(index));
    return *this;
  }

  /// Comparison operators
  /// @{
  /// Equality operator
  constexpr bool operator==(const typed_index& other) const {
    return (id() == other.id()) && (index() == other.index());
  };

  /// Comparison operators
  DETRAY_HOST_DEVICE
  friend constexpr auto operator<=>(const typed_index lhs,
                                    const typed_index rhs) noexcept {
    const auto l{lhs.index()};
    const auto r{rhs.index()};
    if (l < r || (l == r && l < r)) {
      return std::strong_ordering::less;
    }
    if (l > r || (l == r && l > r)) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equivalent;
  }
  /// @}

  /// Arithmetic operators
  /// @{
  DETRAY_HOST_DEVICE
  friend typed_index operator+(const typed_index lhs, const typed_index rhs) {
    return typed_index{}.set_id(lhs.id()).set_index(lhs.index() + rhs.index());
  }

  template <typename idx_t>
    requires(std::integral<idx_t> || std::same_as<idx_t, index_type>)
  DETRAY_HOST_DEVICE friend typed_index operator+(const typed_index lhs,
                                                  const idx_t index) {
    return typed_index{}.set_id(lhs.id()).set_index(lhs.index() + index);
  }

  DETRAY_HOST_DEVICE
  friend typed_index operator-(const typed_index lhs, const typed_index rhs) {
    return typed_index{}.set_id(lhs.id()).set_index(lhs.index() - rhs.index());
  }

  template <typename idx_t>
    requires(std::integral<idx_t> || std::same_as<idx_t, index_type>)
  DETRAY_HOST_DEVICE friend typed_index operator-(const typed_index lhs,
                                                  const idx_t& index) {
    return typed_index{}.set_id(lhs.id()).set_index(lhs.index() - index);
  }

  DETRAY_HOST_DEVICE
  typed_index& operator+=(const typed_index rhs) {
    set_index(this->index() + rhs.index());
    return *this;
  }

  template <typename idx_t>
    requires(std::integral<idx_t> || std::same_as<idx_t, index_type>)
  DETRAY_HOST_DEVICE typed_index& operator+=(const idx_t index) {
    set_index(this->index() + index);
    return *this;
  }

  DETRAY_HOST_DEVICE
  typed_index& operator-=(const typed_index rhs) {
    set_index(this->index() - rhs.index());
    return *this;
  }

  template <typename idx_t>
    requires(std::integral<idx_t> || std::same_as<idx_t, index_type>)
  DETRAY_HOST_DEVICE typed_index& operator-=(const idx_t index) {
    set_index(this->index() - index);
    return *this;
  }
  /// @}

  /// Only make the prefix operator available
  DETRAY_HOST_DEVICE
  typed_index& operator++() {
    if constexpr (std::integral<index_type>) {
      set_index(this->index() + static_cast<index_type>(1));
    } else {
      set_index(this->index() + 1u);
    }

    return *this;
  }

  /// Shift the contained index.
  template <std::integral idx_t>
  DETRAY_HOST_DEVICE constexpr typed_index& shift(idx_t s) {
    if constexpr (std::integral<index_type>) {
      set_index(this->index() + static_cast<index_type>(s));
    } else {
      set_index(this->index().shift(s));
    }

    return *this;
  }

  /// Check whether the link is valid to use.
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid() const noexcept {
    return encoder::template is_invalid<id_mask, index_mask>(m_value);
  }

  /// Check whether the type id is invalid.
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid_id() const noexcept {
    return encoder::template is_invalid<id_mask>(m_value);
  }

  /// Check whether the index is invalid.
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid_index() const noexcept {
    return encoder::template is_invalid<index_mask>(m_value);
  }

 private:
  /// Print the index
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const typed_index ti) {
    os << "[id = " << ti.id() << "(" << static_cast<value_t>(ti.id()) << ") | "
       << "index = " << ti.index();

    if (ti.is_invalid()) {
      os << " (invalid)";
    }
    os << "]";

    return os;
  }

  /// The encoded value. Default: All bits set to 1 (invalid)
  value_t m_value = ~static_cast<value_t>(0);
};

/// Custom get function for the index_range struct - const
template <std::size_t idx, typename index_t, bool contains_size,
          typename value_t, value_t lower_mask, value_t upper_mask>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const index_range<index_t, contains_size, value_t, lower_mask, upper_mask>&
        index) noexcept {
  static_assert(idx <= 1u, "Index must be 0 or 1");

  if constexpr (idx == 0u) {
    return index.lower();
  } else {
    return index.upper();
  }
}

/// Custom get function for the multi_index struct - const
template <std::size_t idx, typename index_type, std::size_t index_size>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const multi_index<index_type, index_size>& index) noexcept {
  return index[idx];
}

/// Custom get function for the multi_index struct.
template <std::size_t idx, typename index_type, std::size_t index_size>
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    multi_index<index_type, index_size>& index) noexcept {
  return index[idx];
}

/// Custom get function for the typed_index struct. Get the type.
template <std::size_t ID, typename id_type, typename index_type>
  requires(ID == 0)
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const typed_index<id_type, index_type>& index) noexcept {
  return index.id();
}

/// Custom get function for the typed_index struct. Get the index.
template <std::size_t ID, typename id_type, typename index_type>
  requires(ID == 1)
DETRAY_HOST_DEVICE constexpr decltype(auto) get(
    const typed_index<id_type, index_type>& index) noexcept {
  return index.index();
}

/// Overload to check for an invalid index range @param ir
template <typename index_t, bool contains_size, typename value_t,
          value_t lower_mask, value_t upper_mask>
DETRAY_HOST_DEVICE constexpr bool is_invalid_value(
    const index_range<index_t, contains_size, value_t, lower_mask, upper_mask>&
        ir) noexcept {
  return ir.is_invalid();
}

/// Overload to check for an invalid multi index @param mi
template <typename index_type, std::size_t index_size>
DETRAY_HOST_DEVICE constexpr bool is_invalid_value(
    const multi_index<index_type, index_size>& mi) noexcept {
  return mi.is_invalid();
}

/// Overload to check for an invalid typed index link @param ti
template <typename id_t, typename index_t, typename value_t, value_t id_mask,
          value_t index_mask>
DETRAY_HOST_DEVICE constexpr bool is_invalid_value(
    const typed_index<id_t, index_t, value_t, id_mask, index_mask>&
        ti) noexcept {
  return ti.is_invalid();
}

}  // namespace detail

}  // namespace detray
