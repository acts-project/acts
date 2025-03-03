// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <stdexcept>

namespace Acts {

class Surface;

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
class GeometryIdentifier {
 public:
  using Value = std::uint64_t;

  /// Construct from an already encoded value.
  constexpr GeometryIdentifier(Value encoded) : m_value(encoded) {}
  /// Construct default GeometryIdentifier with all values set to zero.
  GeometryIdentifier() = default;
  GeometryIdentifier(GeometryIdentifier&&) = default;
  GeometryIdentifier(const GeometryIdentifier&) = default;
  ~GeometryIdentifier() = default;
  GeometryIdentifier& operator=(GeometryIdentifier&&) = default;
  GeometryIdentifier& operator=(const GeometryIdentifier&) = default;

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

  /// Return the approach identifier.
  constexpr Value passive() const { return getBits(kApproachMask); }

  /// Return the sensitive identifier.
  constexpr Value sensitive() const { return getBits(kSensitiveMask); }

  /// Return the extra identifier
  /// Usage can be experiment-specific, like tagging which kind of detector a
  /// surface object corresponds to, or which subsystem it belongs to
  constexpr Value extra() const { return getBits(kExtraMask); }

  /// Set the volume identifier.
  [[deprecated("Use withVolume() instead")]]
  constexpr GeometryIdentifier& setVolume(Value volume) {
    setBits(kVolumeMask, volume);
    return *this;
  }

  /// Set the boundary identifier.
  [[deprecated("Use withBoundary() instead")]]
  constexpr GeometryIdentifier& setBoundary(Value boundary) {
    setBits(kBoundaryMask, boundary);
    return *this;
  }

  /// Set the layer identifier.
  [[deprecated("Use withLayer() instead")]]
  constexpr GeometryIdentifier& setLayer(Value layer) {
    setBits(kLayerMask, layer);
    return *this;
  }

  /// Set the approach identifier.
  [[deprecated("Use withApproach() instead")]]
  constexpr GeometryIdentifier& setApproach(Value approach) {
    setBits(kApproachMask, approach);
    return *this;
  }

  /// Set the approach identifier - shared with Passive
  [[deprecated("Use withPassive() instead")]]
  constexpr GeometryIdentifier& setPassive(Value approach) {
    setBits(kApproachMask, approach);
    return *this;
  }

  /// Set the sensitive identifier.
  [[deprecated("Use withSensitive() instead")]]
  constexpr GeometryIdentifier& setSensitive(Value sensitive) {
    setBits(kSensitiveMask, sensitive);
    return *this;
  }

  /// Set the extra identifier
  [[deprecated("Use withExtra() instead")]]
  constexpr GeometryIdentifier& setExtra(Value extra) {
    return setBits(kExtraMask, extra);
    return *this;
  }

  /// Return a new identifier with the volume set to @p volume
  /// @param volume the new volume identifier
  /// @return a new identifier with the volume set to @p volume
  [[nodiscard]]
  constexpr GeometryIdentifier withVolume(Value volume) const {
    GeometryIdentifier id = *this;
    id.setBits(kVolumeMask, volume);
    return id;
  }

  /// Return a new identifier with the boundary set to @p boundary
  /// @param boundary the new boundary identifier
  /// @return a new identifier with the boundary set to @p boundary
  [[nodiscard]]
  constexpr GeometryIdentifier withBoundary(Value boundary) const {
    GeometryIdentifier id = *this;
    id.setBits(kBoundaryMask, boundary);
    return id;
  }

  /// Return a new identifier with the layer set to @p layer
  /// @param layer the new layer identifier
  /// @return a new identifier with the layer set to @p layer
  [[nodiscard]]
  constexpr GeometryIdentifier withLayer(Value layer) const {
    GeometryIdentifier id = *this;
    id.setBits(kLayerMask, layer);
    return id;
  }

  /// Return a new identifier with the approach set to @p approach
  /// @param approach the new approach identifier
  /// @return a new identifier with the approach set to @p approach
  [[nodiscard]]
  constexpr GeometryIdentifier withApproach(Value approach) const {
    GeometryIdentifier id = *this;
    id.setBits(kApproachMask, approach);
    return id;
  }

  /// Return a new identifier with the passive set to @p passive
  /// @param passive the new passive identifier
  /// @return a new identifier with the passive set to @p passive
  [[nodiscard]]
  constexpr GeometryIdentifier withPassive(Value passive) const {
    GeometryIdentifier id = *this;
    id.setBits(kApproachMask, passive);
    return id;
  }

  /// Return a new identifier with the sensitive set to @p sensitive
  /// @param sensitive the new sensitive identifier
  /// @return a new identifier with the sensitive set to @p sensitive
  [[nodiscard]]
  constexpr GeometryIdentifier withSensitive(Value sensitive) const {
    GeometryIdentifier id = *this;
    id.setBits(kSensitiveMask, sensitive);
    return id;
  }

  /// Return a new identifier with the extra set to @p extra
  /// @param extra the new extra identifier
  /// @return a new identifier with the extra set to @p extra
  [[nodiscard]]
  constexpr GeometryIdentifier withExtra(Value extra) const {
    GeometryIdentifier id = *this;
    id.setBits(kExtraMask, extra);
    return id;
  }

  /// Get the maximum value for the volume identifier.
  /// @return the maximum value for the volume identifier
  static constexpr Value getMaxVolume() { return getMaxValue(kVolumeMask); }

  /// Get the maximum value for the boundary identifier.
  /// @return the maximum value for the boundary identifier
  static constexpr Value getMaxBoundary() { return getMaxValue(kBoundaryMask); }

  /// Get the maximum value for the layer identifier.
  /// @return the maximum value for the layer identifier
  static constexpr Value getMaxLayer() { return getMaxValue(kLayerMask); }

  /// Get the maximum value for the approach identifier.
  /// @return the maximum value for the approach identifier
  static constexpr Value getMaxApproach() { return getMaxValue(kApproachMask); }

  /// Get the maximum value for the sensitive identifier.
  /// @return the maximum value for the sensitive identifier
  static constexpr Value getMaxSensitive() {
    return getMaxValue(kSensitiveMask);
  }

  /// Get the maximum value for the extra identifier.
  /// @return the maximum value for the extra identifier
  static constexpr Value getMaxExtra() { return getMaxValue(kExtraMask); }

 private:
  // clang-format off
  /// (2^8)-1 = 255 volumes
  static constexpr Value kVolumeMask    = 0xff00000000000000;
  /// (2^8)-1 = 255 boundaries
  static constexpr Value kBoundaryMask  = 0x00ff000000000000;
  /// (2^12)-1 = 4095 layers
  static constexpr Value kLayerMask     = 0x0000fff000000000;
  /// (2^8)-1 = 255 approach surfaces
  static constexpr Value kApproachMask  = 0x0000000ff0000000;
  /// (2^20)-1 = 1048575 sensitive surfaces
  static constexpr Value kSensitiveMask = 0x000000000fffff00;
  /// (2^8)-1 = 255 extra values
  static constexpr Value kExtraMask     = 0x00000000000000ff;
  // clang-format on

  Value m_value = 0;

  /// Extract the bit shift necessary to access the masked values.
  static constexpr int extractShift(Value mask) {
    // use compiler builtin to extract the number of trailing bits from the
    // mask. the builtin should be available on all supported compilers.
    // need unsigned long long version (...ll) to ensure std::uint64_t
    // compatibility.
    // WARNING undefined behaviour for mask == 0 which we should not have.
    return __builtin_ctzll(mask);
  }

  constexpr static Value getMaxValue(Value mask) {
    return mask >> extractShift(mask);
  }

  /// Extract the masked bits from the encoded value.
  constexpr Value getBits(Value mask) const {
    return (m_value & mask) >> extractShift(mask);
  }
  /// Set the masked bits to id in the encoded value.
  constexpr GeometryIdentifier& setBits(Value mask, Value id) {
    if (id > getMaxValue(mask)) {
      throw std::invalid_argument(
          "Value " + std::to_string(id) + " exceeds maximum value " +
          std::to_string(getMaxValue(mask)) + " for this field");
    }

    m_value = (m_value & ~mask) | ((id << extractShift(mask)) & mask);
    // return *this here that we need to write fewer lines in the setXXX
    // methods
    return *this;
  }

  friend constexpr bool operator==(GeometryIdentifier lhs,
                                   GeometryIdentifier rhs) {
    return lhs.m_value == rhs.m_value;
  }

  friend constexpr bool operator<(GeometryIdentifier lhs,
                                  GeometryIdentifier rhs) {
    return lhs.m_value < rhs.m_value;
  }
};

std::ostream& operator<<(std::ostream& os, GeometryIdentifier id);

/// Base class for hooks that can be used to modify the Geometry Identifier
/// during construction. Most common use case is setting the extra bit fields.
struct GeometryIdentifierHook {
  virtual ~GeometryIdentifierHook() = default;
  virtual Acts::GeometryIdentifier decorateIdentifier(
      Acts::GeometryIdentifier identifier, const Acts::Surface& surface) const;
};

}  // namespace Acts
// specialize std::hash so GeometryIdentifier can be used e.g. in an
// unordered_map
namespace std {
template <>
struct hash<Acts::GeometryIdentifier> {
  auto operator()(Acts::GeometryIdentifier gid) const noexcept {
    return std::hash<Acts::GeometryIdentifier::Value>()(gid.value());
  }
};
}  // namespace std
