// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryObject.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYOBJECT_H
#define ACTS_GEOMETRYUTILS_GEOMETRYOBJECT_H

#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryID.hpp"

namespace Acts {

/// @class GeometryObject
///
/// Base class to provide GeometryID interface:
/// - simple set and get
///
/// It also provides the binningPosition method for
/// Geometry geometrical object to be binned in BinnedArrays
///
class GeometryObject
{
public:
  /// default constructor
  GeometryObject() : m_geoID(0) {}

  /// constructor from a ready-made value
  ///
  /// @param geoID the geometry identifier of the object
  GeometryObject(const GeometryID& geoID) : m_geoID(geoID) {}

  /// assignment operator
  ///
  /// @param geoID the source geoID
  GeometryObject&
  operator=(const GeometryObject& geoID)
  {
    if (&geoID != this) m_geoID = geoID.m_geoID;
    return *this;
  }

  /// Return the value
  /// @return the geometry id by reference
  const GeometryID&
  geoID() const;

  /// Force a binning position method
  ///
  /// @param bValue is the value in which you want to bin
  ///
  /// @return vector 3D used for the binning schema
  virtual const Vector3D
  binningPosition(BinningValue bValue) const = 0;

  /// Implement the binningValue
  ///
  /// @param bValue is the dobule in which you want to bin
  ///
  /// @return float to be used for the binning schema
  double
  binningPositionValue(BinningValue bValue) const;

  /// Set the value
  ///
  /// @param geoID the geometry identifier to be assigned
  void
  assignGeoID(const GeometryID& geoID);

protected:
  GeometryID m_geoID;
};

inline const GeometryID&
GeometryObject::geoID() const
{
  return m_geoID;
}

inline void
GeometryObject::assignGeoID(const GeometryID& geoID)
{
  m_geoID = geoID;
}

inline double
GeometryObject::binningPositionValue(BinningValue bValue) const
{
  // now switch
  switch (bValue) {
  // case x
  case Acts::binX: {
    return binningPosition(bValue).x();
  } break;
  // case y
  case Acts::binY: {
    return binningPosition(bValue).y();
  } break;
  // case z
  case Acts::binZ: {
    return binningPosition(bValue).z();
  } break;
  // case
  case Acts::binR: {
    return binningPosition(bValue).perp();
  } break;
  // do nothing for the default
  default:
    return 0.;
  }
}
}

#endif
