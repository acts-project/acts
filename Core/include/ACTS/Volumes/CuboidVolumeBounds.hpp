// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CuboidVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_CUBOIDVOLUMESBOUNDS_H
#define ACTS_VOLUMES_CUBOIDVOLUMESBOUNDS_H 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"

namespace Acts {

class RectangleBounds;
class Volume;
class Surface;

/// @class CuboidVolumeBounds
///
/// Bounds for a cubical Volume, the decomposeToSurfaces method creates a
/// vector of 6 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///    - negativeFaceXY [0] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at negative \f$ z \f$
///    - positiveFaceXY [1] : Rectangular Acts::PlaneSurface, parallel to \f$ xy
/// \f$ plane at positive \f$ z \f$
///    - negativeFaceXY [2] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - positiveFaceXY [3] : Rectangular Acts::PlaneSurface, attached to \f$ yz
/// \f$ plane at negative \f$ x \f$
///    - negativeFaceXY [4] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at negative \f$ y \f$
///    - positiveFaceXY [5] : Rectangular Acts::PlaneSurface, parallel to \f$ zx
/// \f$ plane at positive \f$ y \f$
///
///  @image html CuboidVolumeBounds_decomp.gif

class CuboidVolumeBounds : public VolumeBounds
{
public:
  /// @enum BoundValues for readability
  enum BoundValues { bv_halfX = 0, bv_halfY = 1, bv_halfZ = 2, bv_length = 3 };

  /// Default Constructor
  CuboidVolumeBounds();

  /// Constructor - the box boundaries
  /// @param hlengthx is the half length of the cube in x
  /// @param hlengthy is the half length of the cube in y
  /// @param hlengthz is the half length of the cube in z
  CuboidVolumeBounds(double hlenghtx, double hlenghty, double hlengthz);

  /// Copy Constructor
  CuboidVolumeBounds(const CuboidVolumeBounds& bobo);

  /// Destructor
  virtual ~CuboidVolumeBounds();

  /// Assignment operator
  CuboidVolumeBounds&
  operator=(const CuboidVolumeBounds& bobo);

  /// Virtual constructor
  CuboidVolumeBounds*
  clone() const override;

  /// This method checks if position in the 3D volume
  /// frame is inside the cylinder
  bool
  inside(const Vector3D&, double tol = 0.) const override;

  /// Method to decompose the Bounds into boundarySurfaces
  /// @note this method is a pure factory the volume is resposible
  /// for the memory management
  const std::vector<const Surface*>
  decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

  /// This method returns the halflength in local x
  double
  halflengthX() const;

  /// This method returns the halflength in local y
  double
  halflengthY() const;

  /// This method returns the halflength in local z
  double
  halflengthZ() const;

  /// Output Method for std::ostream
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// Templated dumpT method
  template <class T>
  T&
  dumpT(T& dt) const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// parallel to local xy plane
  std::shared_ptr<const RectangleBounds>
  faceXYRectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// parallel to local yz plane
  std::shared_ptr<const RectangleBounds>
  faceYZRectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  // parallel to local zx plane
  std::shared_ptr<const RectangleBounds>
  faceZXRectangleBounds() const;

  /// The bound values
  std::vector<TDD_real_t>                m_valueStore;
  std::shared_ptr<const RectangleBounds> m_xyBounds;
  std::shared_ptr<const RectangleBounds> m_yzBounds;
  std::shared_ptr<const RectangleBounds> m_zxBounds;
};

inline CuboidVolumeBounds*
CuboidVolumeBounds::clone() const
{
  return new CuboidVolumeBounds(*this);
}

inline bool
CuboidVolumeBounds::inside(const Vector3D& pos, double tol) const
{
  return (fabs(pos.x()) <= m_valueStore.at(bv_halfX) + tol
          && fabs(pos.y()) <= m_valueStore.at(bv_halfY) + tol
          && fabs(pos.z()) <= m_valueStore.at(bv_halfZ) + tol);
}

inline double
CuboidVolumeBounds::halflengthX() const
{
  return m_valueStore.at(bv_halfX);
}

inline double
CuboidVolumeBounds::halflengthY() const
{
  return m_valueStore.at(bv_halfY);
}

inline double
CuboidVolumeBounds::halflengthZ() const
{
  return m_valueStore.at(bv_halfZ);
}

template <class T>
T&
CuboidVolumeBounds::dumpT(T& dt) const
{
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(7);
  dt << "Acts::CuboidVolumeBounds: (halfX, halfY, halfZ) = ";
  dt << "(" << m_valueStore.at(bv_halfX) << ", " << m_valueStore.at(bv_halfY)
     << ", " << m_valueStore.at(bv_halfZ) << ")";
  return dt;
}
}

#endif  // ACTS_VOLUMES_CUBOIDVOLUMESBOUNDS_H
