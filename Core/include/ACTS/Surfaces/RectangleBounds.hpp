// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangleBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_RECTANGLEBOUNDS_H
#define ACTS_SURFACES_RECTANGLEBOUNDS_H 1

#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @class RectangleBounds
///
/// Bounds for a rectangular, planar surface.
/// The two local coordinates Acts::eLOC_X, Acts::eLOC_Y are for legacy reasons
/// also called @f$ phi @f$ respectively @f$ eta @f$. The orientation
/// with respect to the local surface framce can be seen in the attached
/// illustration.
///
/// @image html RectangularBounds.gif

class RectangleBounds : public PlanarBounds
{
public:
  /// @enum BoundValues for readability 
  enum BoundValues { bv_halfX = 0, bv_halfY = 1, bv_length = 2 };

  /// Default Constructor - deleted
  RectangleBounds() = delete;

  /// Constructor with halflength in x and y
  /// @param halfX halflength in X
  /// @param halfY halflength in Y
  RectangleBounds(double halfX, double halfY);

  /// Copy constructor
  RectangleBounds(const RectangleBounds& recbo) : PlanarBounds(recbo) {}

  /// Destructor
  virtual ~RectangleBounds();

  /// Assignment Operator
  RectangleBounds&
  operator=(const RectangleBounds& recbo);

  /// Virtual constructor
  virtual RectangleBounds*
  clone() const override;

  /// Return the type of the bounds for persistency 
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Rectangle;
  }

  /// @copydoc SurfaceBounds::inside
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bchk) const override;

  /// @copydoc SurfaceBounds::insideLoc0
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const override;

  /// @copydoc SurfaceBounds::insideLoc1
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const override;

  /// @copydoc SurfaceBounds::minDistance 
  virtual double
  minDistance(const Vector2D& lpos) const override;

  /// Return method for the half length in X
  double
  halflengthX() const;

  /// Return method for the half length in Y
  double
  halflengthY() const;

  /// Return the vertices - or, the points of the extremas 
  virtual const std::vector<Vector2D>
  vertices() const override;

  /// Output Method for std::ostream 
  virtual std::ostream&
  dump(std::ostream& sl) const override;

};

inline RectangleBounds*
RectangleBounds::clone() const
{
  return new RectangleBounds(*this);
}

inline bool
RectangleBounds::inside(const Vector2D& lpos, const BoundaryCheck& bchk) const
{
  if (bchk.bcType == 0)
    return RectangleBounds::insideLoc0(lpos, bchk.toleranceLoc0)
      && RectangleBounds::insideLoc1(lpos, bchk.toleranceLoc1);

  // a fast FALSE
  double max_ell = bchk.lCovariance(0, 0) > bchk.lCovariance(1, 1)
      ? bchk.lCovariance(0, 0)
      : bchk.lCovariance(1, 1);
  double limit = bchk.nSigmas * sqrt(max_ell);
  if (!RectangleBounds::inside(lpos, limit, limit)) return false;
  // a fast TRUE
  double min_ell = bchk.lCovariance(0, 0) < bchk.lCovariance(1, 1)
      ? bchk.lCovariance(0, 0)
      : bchk.lCovariance(1, 1);
  limit = bchk.nSigmas * sqrt(min_ell);
  if (RectangleBounds::inside(lpos, limit, limit)) return true;

  // compute KDOP and axes for surface polygon
  std::vector<KDOP>     elementKDOP(4);
  std::vector<Vector2D> elementP(4);
  float                 theta = (bchk.lCovariance(1, 0) != 0
                 && (bchk.lCovariance(1, 1) - bchk.lCovariance(0, 0)) != 0)
      ? .5
          * bchk.FastArcTan(2 * bchk.lCovariance(1, 0)
                            / (bchk.lCovariance(1, 1) - bchk.lCovariance(0, 0)))
      : 0.;
  sincosCache scResult = bchk.FastSinCos(theta);
  ActsMatrixD<2, 2> rotMatrix;
  rotMatrix << scResult.cosC, scResult.sinC, -scResult.sinC, scResult.cosC;
  // ellipse is always at (0,0), surface is moved to ellipse position and then
  // rotated
  Vector2D p;
  p << m_valueStore.at(RectangleBounds::bv_halfX),
      m_valueStore.at(RectangleBounds::bv_halfY);
  elementP.at(0) = (rotMatrix * (p - lpos));
  p << m_valueStore.at(RectangleBounds::bv_halfX),
      -m_valueStore.at(RectangleBounds::bv_halfY);
  elementP.at(1) = (rotMatrix * (p - lpos));
  p << -m_valueStore.at(RectangleBounds::bv_halfX),
      m_valueStore.at(RectangleBounds::bv_halfY);
  elementP.at(2) = (rotMatrix * (p - lpos));
  p << -m_valueStore.at(RectangleBounds::bv_halfX),
      -m_valueStore.at(RectangleBounds::bv_halfY);
  elementP.at(3)             = (rotMatrix * (p - lpos));
  std::vector<Vector2D> axis = {elementP.at(0) - elementP.at(1),
                                elementP.at(0) - elementP.at(2),
                                elementP.at(0) - elementP.at(3),
                                elementP.at(1) - elementP.at(2)};
  bchk.ComputeKDOP(elementP, axis, elementKDOP);
  // compute KDOP for error ellipse
  std::vector<KDOP> errelipseKDOP(4);
  bchk.ComputeKDOP(bchk.EllipseToPoly(3), axis, errelipseKDOP);
  // check if KDOPs overlap and return result
  return bchk.TestKDOPKDOP(elementKDOP, errelipseKDOP);
}

inline bool
RectangleBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  return (fabs(lpos[Acts::eLOC_X])
          < m_valueStore.at(RectangleBounds::bv_halfX) + tol0);
}

inline bool
RectangleBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  return (fabs(lpos[Acts::eLOC_Y])
          < m_valueStore.at(RectangleBounds::bv_halfY) + tol1);
}

inline double
RectangleBounds::halflengthX() const
{
  return m_valueStore.at(RectangleBounds::bv_halfX);
}

inline double
RectangleBounds::halflengthY() const
{
  return m_valueStore.at(RectangleBounds::bv_halfY);
}

inline const std::vector<Vector2D>
RectangleBounds::vertices() const
{
  // create the return vector
  std::vector<Vector2D> vertices;
  // fill the vertices
  vertices.reserve(4);
  vertices.push_back(
      Vector2D(m_valueStore.at(RectangleBounds::bv_halfX),
               -m_valueStore.at(RectangleBounds::bv_halfY)));  // [0]
  vertices.push_back(
      Vector2D(m_valueStore.at(RectangleBounds::bv_halfX),
               m_valueStore.at(RectangleBounds::bv_halfY)));  // [1]
  vertices.push_back(
      Vector2D(-m_valueStore.at(RectangleBounds::bv_halfX),
               m_valueStore.at(RectangleBounds::bv_halfY)));  // [1]
  vertices.push_back(
      Vector2D(-m_valueStore.at(RectangleBounds::bv_halfX),
               -m_valueStore.at(RectangleBounds::bv_halfY)));  // [3]
  return vertices;
}

}  // end of namespace

#endif  // ACTS_SURFACES_RECTANGLEBOUNDS_H
