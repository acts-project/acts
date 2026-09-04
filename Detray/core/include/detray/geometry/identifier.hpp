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
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/bit_encoder.hpp"

// System include(s)
#include <cstdint>
#include <ostream>
#include <utility>

namespace detray::geometry {

/// @brief Unique identifier for geometry objects, based on the ACTS
/// GeometryIdentifier
///
/// Encodes the volume index, the type of surface (portal, sensitive, passive
/// etc), an index to identify a surface in a geometry accelerator structure,
/// as well as two extra bytes that can be used to tag surfaces arbitrarily.
///
/// @note the detray identifier is not compatible with the ACTS
/// GeometryIdentifier!
///
/// @see
/// https://github.com/acts-project/acts/blob/main/Core/include/Acts/Geometry/GeometryIdentifier.hpp
class identifier {
 public:
  using value_t = std::uint_least64_t;
  using encoder = detail::bit_encoder<value_t>;

  /// Construct from an already encoded value.
  DETRAY_HOST_DEVICE
  constexpr explicit identifier(value_t encoded) : m_value(encoded) {}

  /// Construct default identifier with all values set to 1.
  constexpr identifier() = default;

  /// @returns the encoded value.
  DETRAY_HOST_DEVICE
  constexpr value_t value() const noexcept { return m_value; }

  /// @returns the volume index.
  DETRAY_HOST_DEVICE
  constexpr dindex volume() const {
    return static_cast<dindex>(
        encoder::template get_bits<k_volume_mask>(m_value));
  }

  /// @returns the surface id.
  DETRAY_HOST_DEVICE
  constexpr surface_id id() const {
    return static_cast<surface_id>(
        encoder::template get_bits<k_id_mask>(m_value));
  }

  /// @returns the surface index.
  DETRAY_HOST_DEVICE
  constexpr dindex index() const {
    return static_cast<dindex>(
        encoder::template get_bits<k_index_mask>(m_value));
  }

  /// @returns the transform index
  DETRAY_HOST_DEVICE
  constexpr dindex transform() const {
    return static_cast<dindex>(
        encoder::template get_bits<k_transform_mask>(m_value));
  }

  /// @returns the extra identifier
  DETRAY_HOST_DEVICE
  constexpr dindex extra() const {
    return static_cast<dindex>(
        encoder::template get_bits<k_extra_mask>(m_value));
  }

  /// Set the volume index.
  DETRAY_HOST_DEVICE
  constexpr identifier& set_volume(value_t volume) {
    encoder::template set_bits<k_volume_mask>(m_value, volume);
    return *this;
  }

  /// Set the surface id.
  DETRAY_HOST_DEVICE
  constexpr identifier& set_id(surface_id id) {
    encoder::template set_bits<k_id_mask>(
        m_value, static_cast<value_t>(
                     static_cast<std::underlying_type_t<surface_id>>(id)));
    return *this;
  }

  /// Set the surface index.
  DETRAY_HOST_DEVICE
  constexpr identifier& set_index(value_t index) {
    encoder::template set_bits<k_index_mask>(m_value, index);
    return *this;
  }

  /// Set the transform index.
  DETRAY_HOST_DEVICE
  constexpr identifier& set_transform(value_t index) {
    encoder::template set_bits<k_transform_mask>(m_value, index);
    return *this;
  }

  /// Set the extra identifier.
  DETRAY_HOST_DEVICE
  constexpr identifier& set_extra(value_t extra) {
    encoder::template set_bits<k_extra_mask>(m_value, extra);
    return *this;
  }

  /// Check whether the identifier is valid to use.
  /// @note The extra bits are allowed to be invalid and will not be checked
  DETRAY_HOST_DEVICE
  constexpr bool is_invalid() const noexcept {
    return encoder::template is_invalid<k_volume_mask, k_id_mask, k_index_mask,
                                        k_transform_mask>(m_value);
  }

  bool operator==(const identifier& rhs) const = default;

  /// Comparison operators
  DETRAY_HOST_DEVICE
  friend constexpr auto operator<=>(const identifier lhs,
                                    const identifier rhs) noexcept {
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

 private:
  // clang-format off
    static constexpr value_t k_volume_mask    = 0xfff0000000000000; // (2^12)-1 = 4095 volumes
    static constexpr value_t k_id_mask        = 0x000f000000000000; // (2^4)-1 = 15 surface categories
    static constexpr value_t k_index_mask     = 0x0000fffff8000000; // (2^21)-1 = 2097151 surfaces
    static constexpr value_t k_transform_mask = 0x0000000007ffffc0; // (2^21)-1 = 2097151 transforms
    static constexpr value_t k_extra_mask     = 0x000000000000003f; // (2^6)-1 = 63 extra tags
  // clang-format on

  // The masks together must cover all bits
  static_assert(k_volume_mask + k_id_mask + k_index_mask + k_transform_mask +
                    k_extra_mask ==
                ~static_cast<value_t>(0));

  /// The encoded value. Default: All bits set to 1 (invalid)
  value_t m_value{~static_cast<value_t>(0)};

  /// Print the surface identifier
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os, const identifier c) {
    if (c.is_invalid()) {
      os << "INVALID: ";
    }

    os << "vol = " << c.volume();
    os << " | id = " << c.id() << "(" << static_cast<int>(c.id()) << ")";
    os << " | index = " << c.index();
    os << " | trf = " << c.transform();
    os << " | extra = " << c.extra();

    return os;
  }
};

}  // namespace detray::geometry

namespace std {

// specialize std::hash so identifier can be used e.g. in an unordered_map
template <>
struct hash<detray::geometry::identifier> {
  auto operator()(detray::geometry::identifier gid) const noexcept {
    return std::hash<detray::geometry::identifier::value_t>()(gid.value());
  }
};

}  // namespace std
