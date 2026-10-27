// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <cstdint>
#include <type_traits>

namespace detray::detail {

/// @brief Sets masked bits to a given value.
///
/// Given a mask and an input value, the corresponding bits are set on a target
/// value. Used e.g. in the surface identifier.
///
/// @see
/// https://github.com/acts-project/acts/blob/main/Core/include/Acts/Geometry/GeometryIdentifier.hpp
template <typename value_t = std::uint64_t>
class bit_encoder {
 public:
  bit_encoder() = delete;

  /// Check whether `v` sets every bit in any of the supplied invalid masks.
  ///
  /// @tparam masks masks whose all-ones encoding is reserved as invalid
  /// @param v encoded value to inspect
  /// @return true if at least one mask is fully set in `v`
  template <value_t... masks>
  DETRAY_HOST_DEVICE static constexpr bool is_invalid(value_t v) noexcept {
    // All bits set to one in the range of a given mask defined as invalid
    return (((v & masks) == masks) || ...);
  }

  /// Extract the field selected by `mask` from an encoded value.
  ///
  /// @tparam mask nonzero bit mask selecting the field
  /// @param v encoded value to inspect
  /// @return selected field shifted down to bit zero
  template <value_t mask>
    requires(mask != static_cast<value_t>(0))
  DETRAY_HOST_DEVICE static constexpr value_t get_bits(
      const value_t v) noexcept {
    // Use integral constant to enforce compile time evaluation of shift
    return (v & mask) >> extract_shift(mask);
  }

  /// Replace the field selected by `mask` with `id` in an encoded value.
  ///
  /// @tparam mask nonzero bit mask selecting the field
  /// @param v encoded value to update
  /// @param id field value to encode
  template <value_t mask>
    requires(mask != static_cast<value_t>(0))
  DETRAY_HOST_DEVICE static constexpr void set_bits(value_t& v,
                                                    const value_t id) noexcept {
    // Use integral constant to enforce compile time evaluation of shift
    v = (v & ~mask) | ((id << extract_shift(mask)) & mask);

    // Make sure, the value 'id' can be safely encoded with the bits
    // specified by 'mask' (unless it is an invalid value which has all bits
    // set to 1)
    assert(((id == ~static_cast<value_t>(0)) ||
            (id == static_cast<value_t>(~static_cast<unsigned int>(0))) ||
            (get_bits<mask>(v) == id)) &&
           "Not enough bits in mask to encode value");
  }

 private:
  /// Extract the bit shift necessary to access the masked values.
  ///
  /// @note undefined behaviour for mask == 0 which we should not have.
  DETRAY_HOST_DEVICE
  static constexpr int extract_shift(value_t mask) noexcept {
    return __builtin_ctzll(mask);
  }
};

}  // namespace detray::detail
