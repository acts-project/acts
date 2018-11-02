// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DoubleTrapezoidVolumeBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

namespace Acts {

class Surface;
class RectangleBounds;
class TrapezoidBounds;
class DiamondBounds;

/// @class DoubleTrapezoidVolumeBounds
///
/// Bounds for a double trapezoidal shaped Volume, the decomposeToSurfaces
/// method
/// creates a
/// vector of 8 surfaces:
///
///  BoundarySurfaceFace [index]:
///
///  - negativeFaceXY     [0] : Diamond Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at negative \f$z\f$
///  - positiveFaceXY     [1] : Diamond Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at positive \f$z\f$
///  - trapezoidFaceAlpha1 [2] : Rectangular  Acts::PlaneSurface,
///                              attached to [0] and [1] at negative \f$ x \f$
/// (associated to alpha1)
///  - trapezoidFaceBeta1  [3] : Rectangular  Acts::PlaneSurface,
///                             attached to [0] and [1] at positive \f$ x \f$
/// (associated to beta1)
///  - trapezoidFaceAlpha2 [5] : Rectangular  Acts::PlaneSurface,
///                              attached to [0] and [1] at negative \f$ x \f$
/// (associated to alpha2)
///  - trapezoidFaceBeta2  [6] : Rectangular  Acts::PlaneSurface,
///                              attached to [0] and [1] at positive \f$ x \f$
/// (associated to beta2)
///  - negativeFaceZX     [4] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at negative \f$y\f$
///  - positiveFaceZX     [5] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at positive \f$y\f$
///
///  @image html DoubleTrapezoidVolumeBounds_decomp.gif

class DoubleTrapezoidVolumeBounds : public VolumeBounds
{
public:
  /// @enum BoundValues for readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_medHalfX = 1,
    bv_maxHalfX = 2,
    bv_halfY1   = 3,
    bv_halfY2   = 4,
    bv_halfZ    = 5,
    bv_alpha1   = 6,
    bv_alpha2   = 7,
    bv_length   = 8
  };

  /// Default Constructor
  DoubleTrapezoidVolumeBounds();

  /// Constructor - the double trapezoid boundaries (
  /// symmetric trapezoid/diamond
  ///
  /// @param minhalex half length in x at minimum y
  /// @param medhalex half length in x aty = 0
  /// @param maxhalex half length in x at maximum x
  /// @param haley1 first half length in y (to negative)
  /// @param haley2 second half length in y (to positive)
  /// @param halez half length in z
  DoubleTrapezoidVolumeBounds(double minhalex,
                              double medhalex,
                              double maxhalex,
                              double haley1,
                              double haley2,
                              double halez);

  /// Copy Constructor
  ///
  /// @param trabo is the source bounds
  DoubleTrapezoidVolumeBounds(const DoubleTrapezoidVolumeBounds& trabo);

  /// Destructor
  ~DoubleTrapezoidVolumeBounds() override;

  /// Assignment operator
  ///
  /// @param trabo is the source bounds
  DoubleTrapezoidVolumeBounds&
  operator=(const DoubleTrapezoidVolumeBounds& trabo);

  /// Virtual constructor
  DoubleTrapezoidVolumeBounds*
  clone() const override;

  /// This method checks if position in the 3D volume frame
  /// is inside the cylinder
  ///
  /// @param pos is the global position to be checked for inside
  /// @param tol is the tolerance parametere
  bool
  inside(const Vector3D& pos, double tol = 0.) const override;

  /// decompose into boundary surfaces
  ///
  /// @param transformPtr is the transform applied by the volume
  ///
  /// @return a vector of surfaces to be used as boundary surfaces
  std::vector<std::shared_ptr<const Surface>>
  decomposeToSurfaces(
      std::shared_ptr<const Transform3D> transformPtr) const override;

  /// This method returns the X halflength at minimal Y
  double
  minHalflengthX() const;

  /// This method returns the (maximal) halflength in local x
  double
  medHalflengthX() const;

  /// This method returns the X halflength at maximal Y (local coordinates)
  double
  maxHalflengthX() const;

  /// This method returns the halflength1 in local y
  double
  halflengthY1() const;

  /// This method returns the halflength2 in local y
  double
  halflengthY2() const;

  /// This method returns the halflength in local z
  double
  halflengthZ() const;

  /// This method returns the opening angle in point A (negative local x)
  double
  alpha1() const;

  /// This method returns the opening angle in point A' (negative local x)
  double
  alpha2() const;

  /// Output Method for std::ostream
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// dump method
  ///
  /// @tparam dT is the output stream to be dumped into
  template <class T>
  T&
  dumpT(T& dT) const;

  /// This method returns the associated DoubleTrapezoidBounds of the face
  /// PlaneSurface parallel to local xy plane
  DiamondBounds*
  faceXYDiamondBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// attached to alpha (negative local x)
  RectangleBounds*
  faceAlpha1RectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// attached to alpha (negative local x)
  RectangleBounds*
  faceAlpha2RectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// attached to beta (positive local x)
  RectangleBounds*
  faceBeta1RectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// attached to beta (positive local x)
  RectangleBounds*
  faceBeta2RectangleBounds() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// parallel to local zx plane, negative local y
  RectangleBounds*
  faceZXRectangleBoundsBottom() const;

  /// This method returns the associated RecantleBounds of the face PlaneSurface
  /// parallel to local zx plane, positive local y
  RectangleBounds*
  faceZXRectangleBoundsTop() const;

  std::vector<TDD_real_t> m_valueStore;  ///< the internal store
};

inline DoubleTrapezoidVolumeBounds*
DoubleTrapezoidVolumeBounds::clone() const
{
  return new DoubleTrapezoidVolumeBounds(*this);
}

inline double
DoubleTrapezoidVolumeBounds::minHalflengthX() const
{
  return m_valueStore.at(bv_minHalfX);
}

inline double
DoubleTrapezoidVolumeBounds::medHalflengthX() const
{
  return m_valueStore.at(bv_medHalfX);
}

inline double
DoubleTrapezoidVolumeBounds::maxHalflengthX() const
{
  return m_valueStore.at(bv_maxHalfX);
}

inline double
DoubleTrapezoidVolumeBounds::halflengthY1() const
{
  return m_valueStore.at(bv_halfY1);
}

inline double
DoubleTrapezoidVolumeBounds::halflengthY2() const
{
  return m_valueStore.at(bv_halfY2);
}

inline double
DoubleTrapezoidVolumeBounds::halflengthZ() const
{
  return m_valueStore.at(bv_halfZ);
}

inline double
DoubleTrapezoidVolumeBounds::alpha1() const
{
  return m_valueStore.at(bv_alpha1);
}

inline double
DoubleTrapezoidVolumeBounds::alpha2() const
{
  return m_valueStore.at(bv_alpha2);
}

template <class T>
T&
DoubleTrapezoidVolumeBounds::dumpT(T& dT) const
{
  dT << std::setiosflags(std::ios::fixed);
  dT << std::setprecision(5);
  dT << "Acts::DoubleTrapezoidVolumeBounds: (minhalfX, medhalfX, maxhalfX, "
        "halfY1, halfY2, halfZ) = ";
  dT << "(" << m_valueStore.at(bv_minHalfX) << ", "
     << m_valueStore.at(bv_medHalfX) << ", " << m_valueStore.at(bv_maxHalfX);
  dT << ", " << m_valueStore.at(bv_halfY1) << ", " << m_valueStore.at(bv_halfY2)
     << ", " << m_valueStore.at(bv_halfZ) << ")";
  return dT;
}
}