// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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
  /// constructor from a ready-made value
  GeometryObject(const GeometryID& geoID = GeometryID()) : m_geoID() {}
  /// return the value
  const GeometryID&
  geoID() const;

  /// set the value
  void
  assignGeoID(const GeometryID& geoID) const;

  /// force a binning position method
  virtual const Vector3D
  binningPosition(BinningValue bValue) const = 0;

  /// implement the binningValue
  double
  binningPositionValue(BinningValue bValue) const;

protected:
  mutable GeometryID m_geoID;
};

inline const GeometryID&
GeometryObject::geoID() const
{
  return m_geoID;
}

inline void
GeometryObject::assignGeoID(const GeometryID& geoID) const
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
