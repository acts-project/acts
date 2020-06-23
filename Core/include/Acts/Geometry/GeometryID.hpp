// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <functional>
#include <iosfwd>

namespace Acts {

/// Identifier for geometry nodes within the geometry hierarchy.
///
/// An identifier can be split into the following components. They define
/// a hierarchy of objects starting from the high-level volumes:
///
/// - Volume
/// - Boundary surfaces (for a volume)
/// - Layers (confined within a volume)
/// - Approach surfaces (for a layer)
/// - Sensitive surfaces (confined to a layer, also called modules)
///
class GeometryID {
 public:
  using Value = uint64_t;

  /// Construct from an already encoded value.
  constexpr GeometryID(Value encoded) : m_value(encoded) {}
  /// Construct default GeometryID with all values set to zero.
  GeometryID() = default;
  GeometryID(GeometryID&&) = default;
  GeometryID(const GeometryID&) = default;
  ~GeometryID() = default;
  GeometryID& operator=(GeometryID&&) = default;
  GeometryID& operator=(const GeometryID&) = default;

  /// Return the encoded value.
  constexpr Value value() const { return m_value; }

  /// Return the volume identifier.
  constexpr Value volume() const { return getBits(kVolumeMask); }
  /// Return the boundary identifier.
  constexpr Value boundary() const { return getBits(kBoundaryMask); }
  /// Return the layer identifier.
  constexpr Value layer() const { return getBits(kLayerMask); }
  /// Return the approach identifier.
  constexpr Value approach() const { return getBits(kApproachMask); }
  /// Return the sensitive identifier.
  constexpr Value sensitive() const { return getBits(kSensitiveMask); }

  /// Set the volume identifier.
  constexpr GeometryID& setVolume(Value volume) {
    return setBits(kVolumeMask, volume);
  }
  /// Set the boundary identifier.
  constexpr GeometryID& setBoundary(Value boundary) {
    return setBits(kBoundaryMask, boundary);
  }
  /// Set the layer identifier.
  constexpr GeometryID& setLayer(Value layer) {
    return setBits(kLayerMask, layer);
  }
  /// Set the approach identifier.
  constexpr GeometryID& setApproach(Value approach) {
    return setBits(kApproachMask, approach);
  }
  /// Set the sensitive identifier.
  constexpr GeometryID& setSensitive(Value sensitive) {
    return setBits(kSensitiveMask, sensitive);
  }

 private:
  // (2^8)-1 = 255 volumes
  static constexpr Value kVolumeMask = 0xff00000000000000;
  // (2^8)-1 = 255 boundaries
  static constexpr Value kBoundaryMask = 0x00ff000000000000;
  // (2^12)-1 = 4096 layers
  static constexpr Value kLayerMask = 0x0000fff000000000;
  // (2^8)-1 = 255 approach surfaces
  static constexpr Value kApproachMask = 0x0000000ff0000000;
  // (2^28)-1 sensitive surfaces
  static constexpr Value kSensitiveMask = 0x000000000fffffff;

  Value m_value = 0;

  /// Extract the bit shift necessary to access the masked values.
  static constexpr int extractShift(Value mask) {
    // use compiler builtin to extract the number of trailing bits from the
    // mask. the builtin should be available on all supported compilers.
    // need unsigned long long version (...ll) to ensure uint64_t compatibility.
    // WARNING undefined behaviour for mask == 0 which we should not have.
    return __builtin_ctzll(mask);
  }
  /// Extract the masked bits from the encoded value.
  constexpr Value getBits(Value mask) const {
    return (m_value & mask) >> extractShift(mask);
  }
  /// Set the masked bits to id in the encoded value.
  constexpr GeometryID& setBits(Value mask, Value id) {
    m_value = (m_value & ~mask) | ((id << extractShift(mask)) & mask);
    // return *this here so we need to write less lines in the set... methods
    return *this;
  }

  friend constexpr bool operator==(GeometryID lhs, GeometryID rhs) {
    return lhs.m_value == rhs.m_value;
  }
  friend constexpr bool operator!=(GeometryID lhs, GeometryID rhs) {
    return lhs.m_value != rhs.m_value;
  }
  friend constexpr bool operator<(GeometryID lhs, GeometryID rhs) {
    return lhs.m_value < rhs.m_value;
  }
};

std::ostream& operator<<(std::ostream& os, GeometryID id);

}  // namespace Acts

// specialize std::hash so GeometryId can be used e.g. in an unordered_map
namespace std {
template <>
struct hash<Acts::GeometryID> {
  auto operator()(Acts::GeometryID gid) const noexcept {
    return std::hash<Acts::GeometryID::Value>()(gid.value());
  }
};
}  // namespace std
