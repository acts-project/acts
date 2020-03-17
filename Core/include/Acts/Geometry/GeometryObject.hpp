// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @class GeometryObject
///
/// Base class to provide GeometryID interface:
/// - simple set and get
///
/// It also provides the binningPosition method for
/// Geometry geometrical object to be binned in BinnedArrays
///
class GeometryObject {
 public:
  /// Defaulted construrctor
  GeometryObject() = default;

  /// Defaulted copy constructor
  GeometryObject(const GeometryObject&) = default;

  /// Constructor from a value
  ///
  /// @param geoID the geometry identifier of the object
  GeometryObject(const GeometryID& geoID) : m_geoID(geoID) {}

  /// Assignment operator
  ///
  /// @param geoID the source geoID
  GeometryObject& operator=(const GeometryObject& geoID) {
    if (&geoID != this) {
      m_geoID = geoID.m_geoID;
    }
    return *this;
  }

  /// @return the geometry id by reference
  const GeometryID& geoID() const;

  /// Force a binning position method
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the value in which you want to bin
  ///
  /// @return vector 3D used for the binning schema
  virtual const Vector3D binningPosition(const GeometryContext& gctx,
                                         BinningValue bValue) const = 0;

  /// Implement the binningValue
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the dobule in which you want to bin
  ///
  /// @return float to be used for the binning schema
  virtual double binningPositionValue(const GeometryContext& gctx,
                                      BinningValue bValue) const;

  /// Set the value
  ///
  /// @param geoID the geometry identifier to be assigned
  void assignGeoID(const GeometryID& geoID);

 protected:
  GeometryID m_geoID;
};

inline const GeometryID& GeometryObject::geoID() const {
  return m_geoID;
}

inline void GeometryObject::assignGeoID(const GeometryID& geoID) {
  m_geoID = geoID;
}

inline double GeometryObject::binningPositionValue(const GeometryContext& gctx,
                                                   BinningValue bValue) const {
  return VectorHelpers::cast(binningPosition(gctx, bValue), bValue);
}
}  // namespace Acts
