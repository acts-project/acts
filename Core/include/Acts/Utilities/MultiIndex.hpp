// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cassert>
#include <climits>
#include <ostream>
#include <type_traits>
#include <utility>

namespace Acts {

/// A set of (hierarchical) indices bitpacked into a single value.
///
/// The underlying value is split into blocks of bits with variable size.
/// Each block is a level within the index hierarchy and can be set and
/// retrieved separately. The encoded MultiIndex can be ordered and compared
/// for equality. The ordering follows the hierarchy, i.e. indices are
/// first ordered by the highest level, then within the highest level by the
/// second level and so on.
template <typename T, std::size_t... BitsPerLevel>
class MultiIndex {
 public:
  static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>,
                "The underlying storage type must be an unsigned integer");
  static_assert(0 < sizeof...(BitsPerLevel),
                "At least one level must be defined");
  static_assert((sizeof(T) * CHAR_BIT) == (... + BitsPerLevel),
                "The sum of bits per level must match the underlying storage");

  /// The type of their underlying storage value.
  using Value = T;
  /// Number of levels in the multi-index hierarchy
  static constexpr std::size_t kNumLevels = sizeof...(BitsPerLevel);

  /// Construct a MultiIndex with all levels set to zero.
  /// @return MultiIndex with all levels initialized to zero
  static constexpr MultiIndex Zeros() { return MultiIndex(0u); }
  /// Construct a MultiIndex from values for multiple level.
  ///
  /// This functionality must be implemented as a static, named constructor
  /// to avoid confusion with other constructors. If it would be implemented
  /// as a regular constructor, constructing a MultiIndex from a single
  /// encoded value and encoding only the first level would have the same
  /// signature and could not be distinguished.
  /// @param us Values for each index level to encode
  /// @return MultiIndex encoded with the provided level values
  template <typename... Us>
  static constexpr MultiIndex Encode(Us&&... us) {
    static_assert(sizeof...(Us) <= kNumLevels,
                  "Can only encode as many levels as in the MultiIndex");

    MultiIndex index(0u);
    std::size_t lvl = 0;
    for (Value val : std::array<Value, sizeof...(Us)>{us...}) {
      index.set(lvl++, val);
    }
    return index;
  }

  /// Construct a MultiIndex from an already encoded value.
  /// @param encoded Pre-encoded multi-index value
  explicit constexpr MultiIndex(Value encoded) : m_value(encoded) {}
  /// Construct a default MultiIndex with undefined values for each level.
  MultiIndex() = default;
  /// @brief Copy constructor
  MultiIndex(const MultiIndex&) = default;
  /// @brief Non-const copy constructor
  MultiIndex(MultiIndex&) = default;
  /// @brief Copy assignment operator
  /// @return Reference to this MultiIndex
  MultiIndex& operator=(const MultiIndex&) = default;
  /// @brief Move assignment operator
  /// @return Reference to this MultiIndex
  MultiIndex& operator=(MultiIndex&&) noexcept = default;
  /// Allow setting the MultiIndex from an already encoded value.
  /// @param encoded Pre-encoded multi-index value to assign
  /// @return Reference to this MultiIndex for chaining
  constexpr MultiIndex& operator=(Value encoded) {
    m_value = encoded;
    return *this;
  }

  /// Get the encoded value of all index levels.
  /// @return The complete encoded multi-index value
  constexpr Value value() const { return m_value; }
  /// Get the value for the index level.
  /// @param lvl Level index to retrieve
  /// @return Value stored at the specified level
  constexpr Value level(std::size_t lvl) const {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    return (m_value >> shift(lvl)) & mask(lvl);
  }
  /// Set the value of the index level.
  /// @param lvl Level index to set
  /// @param val Value to set for the specified level
  /// @return Reference to this MultiIndex for chaining
  constexpr MultiIndex& set(std::size_t lvl, Value val) {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    if (val > maxValue(lvl)) {
      throw std::out_of_range(
          "Value " + std::to_string(val) + " for index level " +
          std::to_string(lvl) +
          " exceeds allowed range (max=" + std::to_string(maxValue(lvl)) + ")");
    }
    // mask of valid bits at the encoded positions for the index level
    Value shiftedMask = (mask(lvl) << shift(lvl));
    // value of the index level shifted to its encoded position
    Value shiftedValue = (val << shift(lvl));
    // combine existing values with the value for the index level
    m_value = (m_value & ~shiftedMask) | (shiftedValue & shiftedMask);
    return *this;
  }

  /// @brief Return the maximum allowed value for a specific index level
  /// @param lvl Level index to query maximum value for
  /// @return Maximum value that can be stored at the specified level
  constexpr std::size_t maxValue(std::size_t lvl) const {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    return (1 << s_bits.at(lvl)) - 1;
  }

  /// Create index with the selected level increased and levels below zeroed.
  /// @param lvl Level to increment for creating the sibling
  /// @return New MultiIndex with the next sibling at the specified level
  constexpr MultiIndex makeNextSibling(std::size_t lvl) const {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    // remove lower levels by shifting the upper levels to the left edge
    Value upper = (m_value >> shift(lvl));
    // increase to create sibling and shift back to zero lower levels again
    return MultiIndex{(upper + 1u) << shift(lvl)};
  }
  /// Create index with every level below the selected level maximized.
  /// @param lvl Level below which all levels are maximized
  /// @return New MultiIndex representing the last descendant
  constexpr MultiIndex makeLastDescendant(std::size_t lvl) const {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    // mask everything below the selected level
    Value maskLower = (Value{1u} << shift(lvl)) - 1u;
    // replace the masked lower levels w/ ones
    return MultiIndex{(m_value & ~maskLower) | maskLower};
  }

  /// Get the number of bits for the associated level
  /// @param lvl Level index to query
  /// @return Number of bits allocated for the specified level
  static constexpr std::size_t bits(std::size_t lvl) {
    assert((lvl < kNumLevels) && "Index level outside allowed range");
    return s_bits[lvl];
  }

 private:
  // per-level mask and right-most bit position for shifting
  static constexpr std::array<std::size_t, kNumLevels> s_bits{BitsPerLevel...};
  static constexpr std::size_t shift(std::size_t lvl) {
    std::size_t s = 0u;
    // sum up all bits below the requested level
    for (std::size_t i = (lvl + 1); i < s_bits.size(); ++i) {
      s += s_bits[i];
    }
    return s;
  }
  static constexpr Value mask(std::size_t lvl) {
    return (Value{1u} << s_bits[lvl]) - 1u;
  }

  Value m_value;

  friend constexpr bool operator<(MultiIndex lhs, MultiIndex rhs) {
    return lhs.m_value < rhs.m_value;
  }

  friend constexpr bool operator==(MultiIndex lhs, MultiIndex rhs) {
    return lhs.m_value == rhs.m_value;
  }

  friend inline std::ostream& operator<<(std::ostream& os, MultiIndex idx) {
    // one level is always defined
    os << idx.level(0u);
    for (std::size_t lvl = 1; lvl < kNumLevels; ++lvl) {
      os << '|' << idx.level(lvl);
    }
    return os;
  }
};

}  // namespace Acts

// specialize std::hash so MultiIndex can be used e.g. in an unordered_map
namespace std {
template <typename Storage, std::size_t... BitsPerLevel>
struct hash<Acts::MultiIndex<Storage, BitsPerLevel...>> {
  auto operator()(
      Acts::MultiIndex<Storage, BitsPerLevel...> idx) const noexcept {
    return std::hash<Storage>()(idx.value());
  }
};
}  // namespace std
