// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DoubleTrapezoidVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_DOUBLETRAPEZOIDVOLUMESBOUNDS_H
#define ACTS_VOLUMES_DOUBLETRAPEZOIDVOLUMESBOUNDS_H 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"

namespace Acts {

class Surface;
class RectangleBounds;
class TrapezoidBounds;
class DiamondBounds;

/**
 @class DoubleTrapezoidVolumeBounds

 Bounds for a double trapezoidal shaped Volume, the decomposeToSurfaces method
 creates a
 vector of 8 surfaces:

  BoundarySurfaceFace [index]:

      - negativeFaceXY     [0] : Diamond Acts::PlaneSurface,
                                 parallel to \f$ xy \f$ plane at negative \f$ z
 \f$
      - positiveFaceXY     [1] : Diamond Acts::PlaneSurface,
                                 parallel to \f$ xy \f$ plane at positive \f$ z
 \f$
      - trapezoidFaceAlpha1 [2] : Rectangular  Acts::PlaneSurface,
                                 attached to [0] and [1] at negative \f$ x \f$
 (associated to alpha1)
      - trapezoidFaceBeta1  [3] : Rectangular  Acts::PlaneSurface,
                                 attached to [0] and [1] at positive \f$ x \f$
 (associated to beta1)
      - trapezoidFaceAlpha2 [5] : Rectangular  Acts::PlaneSurface,
                                 attached to [0] and [1] at negative \f$ x \f$
 (associated to alpha2)
      - trapezoidFaceBeta2  [6] : Rectangular  Acts::PlaneSurface,
                                 attached to [0] and [1] at positive \f$ x \f$
 (associated to beta2)
      - negativeFaceZX     [4] : Rectangular  Acts::PlaneSurface,
                                 parallel to \f$ zx \f$ plane at negative \f$ y
 \f$
      - positiveFaceZX     [5] : Rectangular  Acts::PlaneSurface,
                                 parallel to \f$ zx \f$ plane at positive \f$ y
 \f$

  @image html DoubleTrapezoidVolumeBounds_decomp.gif

  */

class DoubleTrapezoidVolumeBounds : public VolumeBounds
{
public:
  /** @enum BoundValues for readability */
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

  /**Default Constructor*/
  DoubleTrapezoidVolumeBounds();

  /**Constructor - the double trapezoid boundaries (symmetric trapezoid/diamond)
   */
  DoubleTrapezoidVolumeBounds(double minhlenghtx,
                              double medhlengthx,
                              double maxhlengthx,
                              double hlenghty1,
                              double hlenghty2,
                              double hlengthz);

  /**Copy Constructor */
  DoubleTrapezoidVolumeBounds(const DoubleTrapezoidVolumeBounds& bobo);

  /**Destructor */
  virtual ~DoubleTrapezoidVolumeBounds();

  /**Assignment operator*/
  DoubleTrapezoidVolumeBounds&
  operator=(const DoubleTrapezoidVolumeBounds& bobo);

  /**Virtual constructor */
  DoubleTrapezoidVolumeBounds*
  clone() const override;

  /**This method checks if position in the 3D volume frame is inside the
   * cylinder*/
  bool
  inside(const Vector3D&, double tol = 0.) const override;

  /** Method to decompose the Bounds into Surfaces */
  const std::vector<const Acts::Surface*>*
  decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

  /**This method returns the X halflength at minimal Y*/
  double
  minHalflengthX() const;

  /**This method returns the (maximal) halflength in local x*/
  double
  medHalflengthX() const;

  /**This method returns the X halflength at maximal Y (local coordinates)*/
  double
  maxHalflengthX() const;

  /**This method returns the halflength1 in local y*/
  double
  halflengthY1() const;

  /**This method returns the halflength2 in local y*/
  double
  halflengthY2() const;

  /**This method returns the halflength in local z*/
  double
  halflengthZ() const;

  /**This method returns the opening angle in point A (negative local x)*/
  double
  alpha1() const;

  /**This method returns the opening angle in point A' (negative local x)*/
  double
  alpha2() const;

  /** Output Method for std::ostream */
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  template <class T>
  T&
  dumpT(T& dT) const;

  /** This method returns the associated DoubleTrapezoidBounds of the face
   * PlaneSurface parallel to local xy plane */
  DiamondBounds*
  faceXYDiamondBounds() const;

  /** This method returns the associated RecantleBounds of the face PlaneSurface
   * attached to alpha (negative local x)*/
  RectangleBounds*
  faceAlpha1RectangleBounds() const;
  RectangleBounds*
  faceAlpha2RectangleBounds() const;

  /** This method returns the associated RecantleBounds of the face PlaneSurface
   * attached to beta (positive local x)*/
  RectangleBounds*
  faceBeta1RectangleBounds() const;
  RectangleBounds*
  faceBeta2RectangleBounds() const;

  /** This method returns the associated RecantleBounds of the face PlaneSurface
   * parallel to local zx plane, negative local y */
  RectangleBounds*
  faceZXRectangleBoundsBottom() const;

  /** This method returns the associated RecantleBounds of the face PlaneSurface
   * parallel to local zx plane, positive local y */
  RectangleBounds*
  faceZXRectangleBoundsTop() const;

  std::vector<TDD_real_t> m_valueStore;
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
  dT << std::setprecision(7);
  dT << "Acts::DoubleTrapezoidVolumeBounds: (minhalfX, medhalfX, maxhalfX, "
        "halfY1, halfY2, halfZ) = ";
  dT << "(" << m_valueStore.at(bv_minHalfX) << ", "
     << m_valueStore.at(bv_medHalfX) << ", " << m_valueStore.at(bv_maxHalfX);
  dT << ", " << m_valueStore.at(bv_halfY1) << ", "
     << m_valueStore.at(bv_halfY2) << ", " << m_valueStore.at(bv_halfZ)
     << ")";
  return dT;
}
}

#endif  // ACTS_VOLUMES_DOUBLETRAPEZOIDVOLUMESBOUNDS_H
