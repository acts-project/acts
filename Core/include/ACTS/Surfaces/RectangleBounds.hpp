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
  ///
  /// @param halfX halflength in X
  /// @param halfY halflength in Y
  RectangleBounds(double halfX, double halfY);

  /// Copy constructor
  ///
  /// @param recbo are the source bounds
  RectangleBounds(const RectangleBounds& recbo) : PlanarBounds(recbo) {}

  /// Destructor
  virtual ~RectangleBounds();

  /// Assignment Operator
  ///
  /// @param recbo are the source bounds
  RectangleBounds&
  operator=(const RectangleBounds& recbo);

  /// Virtual constructor
  virtual RectangleBounds*
  clone() const final override;

  /// Return the type of the bounds for persistency
  virtual BoundsType
  type() const final override
  {
    return SurfaceBounds::Rectangle;
  }

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Inside check for the bounds object with tolerance
  /// checks for first coordinate only.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absolute tolerance parameter
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const final override;

  /// Inside check for the bounds object with tolerance
  /// checks for second coordinate only.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol1 absulote tolerance parameter
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const final override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const final override;

  /// Return method for the half length in X
  double
  halflengthX() const;

  /// Return method for the half length in Y
  double
  halflengthY() const;

  /// Return the vertices - or, the points of the extremas
  virtual const std::vector<Vector2D>
  vertices() const final override;

  // Bounding box representation
  virtual const RectangleBounds&
  boundingBox() const final override;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream for the dump
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

private:
  /// Private helper method
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absulote tolerance parameter on the first coordinate
  /// @param tol1 absulote tolerance parameter on the second coordinate
  ///
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D& lpos, double tol0 = 0., double tol1 = 0.) const;
};

inline RectangleBounds*
RectangleBounds::clone() const
{
  return new RectangleBounds(*this);
}

inline bool
RectangleBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  if (bcheck.bcType == 0)
    return RectangleBounds::insideLoc0(lpos, bcheck.toleranceLoc0)
        && RectangleBounds::insideLoc1(lpos, bcheck.toleranceLoc1);

  // a fast FALSE
  double max_ell = (*bcheck.lCovariance)(0, 0) > (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  double limit = bcheck.nSigmas * sqrt(max_ell);
  if (!RectangleBounds::inside(lpos, limit, limit)) return false;
  // a fast TRUE
  double min_ell = (*bcheck.lCovariance)(0, 0) < (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  limit = bcheck.nSigmas * sqrt(min_ell);
  if (RectangleBounds::inside(lpos, limit, limit)) return true;

  // compute KDOP and axes for surface polygon
  std::vector<KDOP>     elementKDOP(4);
  std::vector<Vector2D> elementP(4);
  float                 theta
      = ((*bcheck.lCovariance)(1, 0) != 0
         && ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)) != 0)
      ? .5
          * bcheck.FastArcTan(
                2 * (*bcheck.lCovariance)(1, 0)
                / ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)))
      : 0.;
  sincosCache scResult = bcheck.FastSinCos(theta);
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
  bcheck.ComputeKDOP(elementP, axis, elementKDOP);
  // compute KDOP for error ellipse
  std::vector<KDOP> errelipseKDOP(4);
  bcheck.ComputeKDOP(bcheck.EllipseToPoly(3), axis, errelipseKDOP);
  // check if KDOPs overlap and return result
  return bcheck.TestKDOPKDOP(elementKDOP, errelipseKDOP);
}

inline bool
RectangleBounds::inside(const Vector2D& lpos, double tol0, double tol1) const
{
  return insideLoc0(lpos, tol0) && insideLoc1(lpos, tol1);
}

inline bool
RectangleBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  return (std::abs(lpos[Acts::eLOC_X])
          < m_valueStore.at(RectangleBounds::bv_halfX) + tol0);
}

inline bool
RectangleBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  return (std::abs(lpos[Acts::eLOC_Y])
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

inline const RectangleBounds&
RectangleBounds::boundingBox() const
{
  return (*this);
}

}  // end of namespace

#endif  // ACTS_SURFACES_RECTANGLEBOUNDS_H
