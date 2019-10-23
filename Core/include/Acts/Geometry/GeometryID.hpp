// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iostream>

#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

using geo_id_value = uint64_t;

/// @class GeometryID
///
///  Identifier for Geometry nodes - packing the
///  - (Sensitive) Surfaces    - uses counting through sensitive surfaces
///  - (Approach)  Surfaces    - uses counting approach surfaces
///  - (Layer)     Surfaces    - uses counting confined layers
///  - (Boundary)  Surfaces    - uses counting through boundary surfaces
///  - Volumes                 - uses counting given by TrackingGeometry

class GeometryID {
 public:
  // clang-format off
  constexpr static geo_id_value volume_mask    = 0xff00000000000000; // 255 volumes
  constexpr static geo_id_value boundary_mask  = 0x00ff000000000000; // 255 boundaries
  constexpr static geo_id_value layer_mask     = 0x0000fff000000000; // 4095 layers
  constexpr static geo_id_value approach_mask  = 0x0000000ff0000000; // 255 approach surfaces
  constexpr static geo_id_value sensitive_mask = 0x000000000fffffff; // (2^28)-1 sensitive surfaces
  // clang-format on

  /// Construct default GeometryID with all values set to zero.
  GeometryID() : m_value(0) {}

  /// constructor from a ready-made value
  ///
  /// @param id_value is the full decoded value of the identifier
  GeometryID(geo_id_value id_value) : m_value(id_value) {}

  // constructor from a shift and a value
  ///
  /// @param id type_id numbered object
  /// @param type_mask is necessary for the decoding
  GeometryID(geo_id_value type_id, geo_id_value type_mask)
      : m_value(ACTS_BIT_ENCODE(type_id, type_mask)) {}

  /// Copy constructor
  ///
  /// @param tddID is the geometry ID that will be copied
  GeometryID(const GeometryID& tddID) = default;

  /// Assignement operator
  ///
  /// @param tddID is the geometry ID that will be assigned
  GeometryID& operator=(const GeometryID& tddID) {
    if (&tddID != this) {
      m_value = tddID.m_value;
    }
    return (*this);
  }

  /// Equality operator
  ///
  /// @param tddID is the geometry ID that will be compared on equality
  bool operator==(const GeometryID& tddID) const {
    return (m_value == tddID.value());
  }

  /// Non-equality operator
  ///
  /// @param tddID is the geometry ID that will be compared on equality
  bool operator!=(const GeometryID& tddID) const { return !operator==(tddID); }

  /// Return the encoded value.
  constexpr geo_id_value value() const { return m_value; }

  /// Return the volume identifier.
  constexpr geo_id_value volume() const {
    return ACTS_BIT_DECODE(m_value, volume_mask);
  }
  /// Return the boundary identifier.
  constexpr geo_id_value boundary() const {
    return ACTS_BIT_DECODE(m_value, boundary_mask);
  }
  /// Return the layer identifier.
  constexpr geo_id_value layer() const {
    return ACTS_BIT_DECODE(m_value, layer_mask);
  }
  /// Return the approach identifier.
  constexpr geo_id_value approach() const {
    return ACTS_BIT_DECODE(m_value, approach_mask);
  }
  /// Return the sensitive identifier.
  constexpr geo_id_value sensitive() const {
    return ACTS_BIT_DECODE(m_value, sensitive_mask);
  }

  /// Set the volume identifier.
  constexpr GeometryID& setVolume(geo_id_value volume) {
    return setBits(volume_mask, volume);
  }
  /// Set the boundary identifier.
  constexpr GeometryID& setBoundary(geo_id_value boundary) {
    return setBits(boundary_mask, boundary);
  }
  /// Set the layer identifier.
  constexpr GeometryID& setLayer(geo_id_value layer) {
    return setBits(layer_mask, layer);
  }
  /// Set the approach identifier.
  constexpr GeometryID& setApproach(geo_id_value approach) {
    return setBits(approach_mask, approach);
  }
  /// Set the sensitive identifier.
  constexpr GeometryID& setSensitive(geo_id_value sensitive) {
    return setBits(sensitive_mask, sensitive);
  }

  /// return the split value as string for debugging
  std::string toString() const;

 private:
  geo_id_value m_value = 0;

  /// Set the subset of bits indicated by the mask
  constexpr GeometryID& setBits(geo_id_value mask, geo_id_value id) {
    m_value = (m_value & ~mask) | ACTS_BIT_ENCODE(id, mask);
    return *this;
  }
};

inline std::string GeometryID::toString() const {
  std::stringstream dstream;
  dstream << "[ " << std::setw(3) << volume();
  dstream << " | " << std::setw(3) << boundary();
  dstream << " | " << std::setw(3) << layer();
  dstream << " | " << std::setw(3) << approach();
  dstream << " | " << std::setw(4) << sensitive() << " ]";
  return dstream.str();
}

/// Overload of operator< | <= | > | >=  for the usage in a map
bool operator<(const GeometryID& one, const GeometryID& two);
bool operator<=(const GeometryID& one, const GeometryID& two);
bool operator>(const GeometryID& one, const GeometryID& two);
bool operator>=(const GeometryID& one, const GeometryID& two);

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const GeometryID& tid);
}  // namespace Acts
