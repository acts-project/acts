// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryID.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYID_H
#define ACTS_GEOMETRYUTILS_GEOMETRYID_H 1

#include <iostream>

namespace Acts {

typedef uint64_t geo_id_value;

/// @class GeometryID
///
///  Identifier for Geometry nodes - packing the
///  - (Sensitive) Surfaces    - uses counting through sensitive surfaces
///  - (Approach)  Surfaces    - uses counting approach surfaces
///  - (Layer)     Surfaces    - uses counting confined layers
///  - (Boundary)  Surfaces    - uses counting through boundary surfaces
///  - Volumes                 - uses counting given by TrackingGeometry

class GeometryID
{
public:
  const static geo_id_value volume_mask     = 0xff00000000000000;
  const static geo_id_value volume_shift    = 56;
  const static geo_id_value volume_range    = 64 - volume_shift;
  const static geo_id_value boundary_mask   = 0x00ff000000000000;
  const static geo_id_value boundary_shift  = 48;
  const static geo_id_value boundary_range  = volume_shift - boundary_shift;
  const static geo_id_value layer_mask      = 0x0000ff0000000000;
  const static geo_id_value layer_shift     = 40;
  const static geo_id_value layer_range     = boundary_shift - layer_shift;
  const static geo_id_value approach_mask   = 0x000000ff00000000;
  const static geo_id_value approach_shift  = 32;
  const static geo_id_value approach_range  = layer_shift - approach_shift;
  const static geo_id_value sensitive_mask  = 0x00000000ffff0000;
  const static geo_id_value sensitive_shift = 16;
  const static geo_id_value sensitive_range = approach_shift - sensitive_shift;
  const static geo_id_value channel_mask    = 0x000000000000ffff;
  const static geo_id_value channel_shift   = 0;
  const static geo_id_value channel_range   = sensitive_shift;

  /// default constructor
  ///
  GeometryID() : m_value(0) {}
  
  /// constructor from a ready-made value
  ///
  /// @param id_value is the full decoded value of the identifier
  GeometryID(geo_id_value id_value) : m_value(id_value) {}
  
  // constructor from a shift and a value
  ///
  /// @param id is numbered object
  /// @param type_shift is the shift necessary for the object type
  GeometryID(geo_id_value id, geo_id_value type_shift)
    : m_value((id << type_shift))
  {}

  /// Copy constructor
  ///
  /// @param tddID is the geometry ID that will be copied
  GeometryID(const GeometryID& tddID) : m_value(tddID.m_value) {}
  
  /// Assignement operator
  ///
  /// @param tddID is the geometry ID that will be assigned
  GeometryID&
  operator=(const GeometryID& tddID)
  {
    if (&tddID != this) m_value = tddID.m_value;
    return (*this);
  }

  /// Add some stuff - a new GeometryID
  ///
  /// @param tddID is the geometry ID that will be added
  GeometryID&
  operator+=(const GeometryID& tddID)
  {
    m_value += tddID.value();
    return (*this);
  }

  /// Add some stuff - a decoded value
  ///
  /// @param add_value is the fully decoded value to be added
  GeometryID&
  operator+=(geo_id_value add_value)
  {
    m_value += add_value;
    return (*this);
  }

  /// return the value
  ///
  /// @param mask is the mask to be applied
  /// @param shift is the according shift to be applied
  geo_id_value
  value(geo_id_value mask = 0, geo_id_value shift = 0) const;

private:
  geo_id_value m_value;
};

inline geo_id_value
GeometryID::value(geo_id_value mask, geo_id_value shift) const
{
  if (mask) return ((m_value & mask) >> shift);
  return m_value;
}

/// Overload of operator< | <= | > | >=  for the usage in a map
bool
operator<(const GeometryID& one, const GeometryID& two);
bool
operator<=(const GeometryID& one, const GeometryID& two);
bool
operator>(const GeometryID& one, const GeometryID& two);
bool
operator>=(const GeometryID& one, const GeometryID& two);

/// Overload of << operator for std::ostream for debug output
std::ostream&
operator<<(std::ostream& sl, const GeometryID& tddID);
}

#endif  // ACTS_GEOMETRYUTILS_GEOMETRYID_H
