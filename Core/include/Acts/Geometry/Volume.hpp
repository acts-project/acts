// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <iosfwd>
#include <memory>

namespace Acts {

class VolumeBounds;

/// @class Volume
///
/// It inhertis of GeometryObject for TDD identification
///
/// Base class for all volumes inside the tracking realm, it defines
/// the interface for inherited Volume classes
/// regarding the geometrical information.

class Volume : public virtual GeometryObject {
 public:
  using BoundingBox = AxisAlignedBoundingBox<Volume, ActsScalar, 3>;

  /// Explicit constructor with shared arguments
  ///
  /// @param transform is the transform to position the volume in 3D space
  /// @param volbounds is the volume boundary definitions
  Volume(const Transform3& transform, std::shared_ptr<VolumeBounds> volbounds);

  /// Copy Constructor - with optional shift
  ///
  /// @param vol is the source volume for the copy
  /// @param shift is the optional shift applied as : shift * vol.transform()
  Volume(const Volume& vol, const Transform3& shift = Transform3::Identity());

  Volume() = delete;
  virtual ~Volume() = default;

  /// Assignment operator
  ///
  /// @param vol is the source volume to be copied
  Volume& operator=(const Volume& vol);

  /// Return methods for geometry transform
  const Transform3& transform() const;

  /// Returns the inverted transform of this volume.
  const Transform3& itransform() const;

  /// returns the center of the volume
  const Vector3& center() const;

  /// Returns const reference to the volume bounds
  const VolumeBounds& volumeBounds() const;

  /// Returns reference to the volume bounds
  /// @warning This is INCOMPATIBLE with Gen1 geometry building:
  ///          changing volumes bounds will invalidate boundary surfaces,
  ///          and the change is not propagated!
  VolumeBounds& volumeBounds();

  /// Set volume bounds and update volume bounding boxes implicitly
  void assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds);

  /// Construct bounding box for this shape
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @return Constructed bounding box pointing to this volume
  BoundingBox boundingBox(const Vector3& envelope = {0, 0, 0}) const;

  /// Construct oriented bounding box for this shape
  /// @note This will build an oriented bounding box with an
  ///       envelope value of (0.05, 0.05, 0.05)mm
  /// @return Constructed oriented bounding box pointing to this volume
  BoundingBox orientedBoundingBox() const;

  /// Inside() method for checks
  ///
  /// @param gpos is the position to be checked
  /// @param tol is the tolerance parameter
  ///
  /// @return boolean indicator if the position is inside
  bool inside(const Vector3& gpos, ActsScalar tol = 0.) const;

  /// The binning position method
  /// - as default the center is given, but may be overloaded
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the binning value schema
  ///
  /// @return vector 3D that can be used for the binning
  Vector3 binningPosition(const GeometryContext& gctx,
                          BinningValue bValue) const override;

 protected:
  Transform3 m_transform;
  Transform3 m_itransform;
  Vector3 m_center;
  std::shared_ptr<VolumeBounds> m_volumeBounds;
};

/**Overload of << operator for std::ostream for debug output*/
std::ostream& operator<<(std::ostream& sl, const Volume& vol);

}  // namespace Acts
