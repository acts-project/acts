// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AbstractVolume.h, ACTS project
///////////////////////////////////////////////////////////////////
#ifndef ACTS_VOLUMES_VOLUME_H
#define ACTS_VOLUMES_VOLUME_H 1

#include <memory>
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
#include "ACTS/Utilities/GeometryStatics.hpp"

namespace Acts {

class VolumeBounds;
typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

/// @class Volume
///
/// It inhertis of GeometryObject for TDD identification
///
/// Base class for all volumes inside the tracking realm, it defines
/// the interface for inherited Volume classes
/// regarding the geometrical information.

class Volume : public virtual GeometryObject
{
public:
  ///  Default constructor
  Volume();

  /// Explicit constructor with shared arguments
  ///
  /// @param htrans is the transform to position the volume in 3D space
  /// @param volBounds is the volume boundary definitions
  Volume(std::shared_ptr<Transform3D> htrans, VolumeBoundsPtr volBounds);

  /// Copy Constructor - with optional shift
  ///
  /// @param vol is the source volume for the copy
  /// @param shift is the optional shift applied after copying
  Volume(const Volume& vol, const Transform3D* shift = nullptr);

  /// Destructor
  virtual ~Volume();

  /// Assignment operator
  ///
  /// @param vol is the source volume to be copied
  Volume&
  operator=(const Volume& vol);

  /// Pseudo-constructor
  virtual Volume*
  clone() const;

  //// Return methods for geometry transform
  const Transform3D&
  transform() const;

  /// returns the center of the volume
  const Vector3D&
  center() const;

  /// returns the volumeBounds()
  const VolumeBounds&
  volumeBounds() const;

  /// Inside() method for checks
  ///
  /// @param gpos is the position to be checked
  /// @param tol is the tolerance parameter
  ///
  /// @return boolean indicator if the position is inside
  bool
  inside(const Vector3D& gpos, double tol = 0.) const;

  /// The binning position method
  /// - as default the center is given, but may be overloaded
  ///
  /// @param bValue is the binning value schema
  ///
  /// @return vector 3D that can be used for the binning
  virtual const Vector3D
  binningPosition(BinningValue bValue) const override;

protected:
  std::shared_ptr<Transform3D> m_transform;
  Vector3D                     m_center;
  VolumeBoundsPtr              m_volumeBounds;
};

inline const Transform3D&
Volume::transform() const
{
  if (m_transform) return (*(m_transform.get()));
  return Acts::s_idTransform;
}

inline const Vector3D&
Volume::center() const
{
  return m_center;
}

inline const VolumeBounds&
Volume::volumeBounds() const
{
  return (*(m_volumeBounds.get()));
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream&
operator<<(std::ostream& sl, const Volume& vol);

}  // end of namespace Acts

#endif  // ACTS_VOLUMES_VOLUME_H
