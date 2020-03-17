// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iomanip>
#include <iostream>
#include <memory>
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"

#ifndef VOLUMEBOUNDS_VALUESTORE_FILL
#define VOLUMEBOUNDS_VALUESTORE_FILL(val) m_valueStore[bv_##val] = val
#endif

#ifndef VOLUMEBOUNDS_VALUESTORE_ACCESS
#define VOLUMEBOUNDS_VALUESTORE_ACCESS(val) \
  double val() const { return m_valueStore[bv_##val]; }
#endif

#ifndef VOLUMEBOUNDS_DERIVED_ACCESS
#define VOLUMEBOUNDS_DERIVED_ACCESS(derived) \
  double derived() const { return m_##derived; }
#endif

namespace Acts {

class Surface;
class Volume;

class VolumeBounds;
using VolumeBoundsPtr = std::shared_ptr<const VolumeBounds>;

using SurfacePtr = std::shared_ptr<const Surface>;
using SurfacePtrVector = std::vector<SurfacePtr>;
/// @class VolumeBounds
///
/// Pure Absract Base Class for Volume bounds.
///
/// Acts::VolumeBounds are a set of up to six confining Surfaces that are stored
/// in a std::vector.
/// Each type of Acts::VolumeBounds has to implement a decomposeToSurfaces() and
/// a inside() method.
///
/// The orientation of the Surfaces are in a way that the normal vector points
/// to the outside world.
///
/// The Volume, retrieving a set of Surfaces from the VolumeBounds, can turn the
/// Surfaces into BoundarySurfaces.
class VolumeBounds {
 public:
  /// Default Constructor*/
  VolumeBounds() = default;

  /// Destructor
  virtual ~VolumeBounds() = default;

  ///  clone() method to make deep copy in Volume copy constructor and for
  /// assigment operator  of the Surface class.
  virtual VolumeBounds* clone() const = 0;

  /// Checking if position given in volume frame is inside
  ///
  /// @param gpos is the global position to be checked
  /// @param tol is the tolerance applied for the inside check
  ///
  /// @return boolean indicating if the position is inside
  virtual bool inside(const Vector3D& gpos, double tol = 0.) const = 0;

  /// Method to decompose the Bounds into Surfaces
  /// the Volume can turn them into BoundarySurfaces
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  ///
  /// @return a vector of surfaces bounding this volume
  virtual SurfacePtrVector decomposeToSurfaces(
      const Transform3D* transform) const = 0;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  virtual Volume::BoundingBox boundingBox(
      const Transform3D* trf = nullptr, const Vector3D& envelope = {0, 0, 0},
      const Volume* entity = nullptr) const = 0;

  /// Binning offset - overloaded for some R-binning types
  ///
  /// @param bValue is the binning schema used
  ///
  /// @return vector 3D to be used for the binning
  virtual Vector3D binningOffset(BinningValue bValue) const;

  /// Binning borders in double
  ///
  /// @param bValue is the binning schema used
  ///
  /// @return float offset to be used for the binning
  virtual double binningBorder(BinningValue bValue) const;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the output stream to be dumped into
  virtual std::ostream& toStream(std::ostream& sl) const = 0;
};

/// Binning offset - overloaded for some R-binning types
inline Vector3D VolumeBounds::binningOffset(
    BinningValue /*bValue*/) const {  // standard offset is 0.,0.,0.
  return Vector3D(0., 0., 0.);
}

inline double VolumeBounds::binningBorder(BinningValue /*bValue*/) const {
  return 0.;
}

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const VolumeBounds& vb);

}  // namespace Acts
