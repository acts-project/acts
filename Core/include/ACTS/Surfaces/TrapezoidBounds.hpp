// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_TRAPEZOIDBOUNDS_H
#define ACTS_SURFACES_TRAPEZOIDBOUNDS_H 1

#include <math.h>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class TrapezoidBounds
///
/// Bounds for a trapezoidal, planar Surface.
///
/// @image html TrapezoidalBounds.gif
///
/// @todo can be speed optimized by calculating kappa/delta and caching it

class TrapezoidBounds : public PlanarBounds
{
public:
  /// @enum BoundValues - for readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_maxHalfX = 1,
    bv_halfY    = 2,
    bv_length   = 3
  };

  /// Trapezoid bounds default constructor is forbidden
  TrapezoidBounds() = delete;

  /// Constructor for symmetric Trapezoid
  ///
  /// @param minhalex minimal half lenght X, definition at negative halflength Y
  /// @param maxhalex maximal half length X, definition at maximum halflength Y
  /// @param haley maximal half length Y - defined at x=0
  TrapezoidBounds(double minhalex, double maxhalex, double haley);

  /// Constructor for arbitrary Trapezoid
  ///
  /// @param minhalex minimal half lenght X, definition at negative halflength Y
  /// @param maxhalex maximal half length X, definition at maximum halflength Y
  /// @param alpha opening angle at @todo check
  /// @param beta opentin angle at @todo check
  TrapezoidBounds(double minhalex, double haley, double alpha, double beta);

  /// Copy constructor
  ///
  /// @param trabo are the source bounds for assignment
  TrapezoidBounds(const TrapezoidBounds& trabo) : PlanarBounds(trabo) {}
  /// Destructor
  virtual ~TrapezoidBounds();

  /// Virtual constructor
  virtual TrapezoidBounds*
  clone() const override;

  /// Return the type of the bounds for persistency
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Trapezoid;
  }

  /// Assignment operator
  TrapezoidBounds&
  operator=(const TrapezoidBounds& sbo);

  ///  This method returns the minimal halflength in X
  /// (first coordinate of local surface frame)
  double
  minHalflengthX() const;

  /// This method returns the maximal halflength in X
  /// (first coordinate of local surface frame)
  double
  maxHalflengthX() const;

  /// This method returns the halflength in Y
  /// (second coordinate of local surface frame)
  double
  halflengthY() const;

  /// This method returns the opening angle alpha in point A
  // (negative local phi)
  double
  alpha() const;

  /// This method returns the opening angle beta in point B
  /// (positive local phi)
  double
  beta() const;

  /// The orientation of the Trapezoid is according to the figure above,
  /// in words: the shorter of the two parallel sides of the trapezoid
  /// intersects
  /// with the negative @f$ y @f$ - axis of the local frame.
  ///
  /// @param lpos is the local position to be checked (carthesian local frame)
  /// @param bcheck is the boundary check directive
  ///
  /// <br>
  /// The cases are:<br>
  /// (0) @f$ y @f$ or @f$ x @f$ bounds are 0 || 0<br>
  /// (1) the local position is outside @f$ y @f$ bounds <br>
  /// (2) the local position is inside @f$ y @f$ bounds, but outside maximum @f$
  /// x
  /// @f$ bounds  <br>
  /// (3) the local position is inside @f$ y @f$ bounds AND inside minimum @f$ x
  /// @f$ bounds <br>
  /// (4) the local position is inside @f$ y @f$ bounds AND inside maximum @f$ x
  /// @f$ bounds, so that
  /// it depends on the @f$ eta @f$ coordinate
  /// (5) the local position fails test of (4) <br>
  ///
  /// The inside check is done using single equations of straight lines and one
  /// has
  /// to take care if a point
  /// lies on the positive @f$ x @f$ half area(I) or the negative one(II).
  /// Denoting
  /// @f$ |x_{min}| @f$ and
  /// @f$ | x_{max} | @f$ as \c minHalfX respectively \c maxHalfX, such as @f$ |
  /// y_{H} | @f$ as \c halfY,
  /// the equations for the straing lines in (I) and (II) can be written as:<br>
  ///  <br>
  /// - (I):  @f$ y = \kappa_{I} x + \delta_{I} @f$ <br>
  /// - (II): @f$ y = \kappa_{II} x + \delta_{II} @f$ ,<br>
  ///  <br>
  /// where @f$  \kappa_{I} = - \kappa_{II} = 2 \frac{y_{H}}{x_{max} - x_{min}}
  /// @f$
  /// <br>
  /// and   @f$  \delta_{I} = \delta_{II} = - \frac{1}{2}\kappa_{I}(x_{max} +
  /// x_{min}) @f$
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const override;

  /// This method checks inside bounds in loc0
  /// @note loc0/loc1 correspond to the natural coordinates of the surface
  /// @note As loc0/loc1 are correlated the single check doesn't make sense :
  ///       check is done on enclosing Rectangle !
  ///
  /// @param lpos is the local position to be checked
  /// @param tol0 is the tolerance applied
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const override;

  /// This method checks inside bounds in loc0
  /// @note loc0/loc1 correspond to the natural coordinates of the surface
  /// @note As loc0/loc1 are correlated the single check doesn't make sense :
  ///   -> check is done on enclosing Rectangle !
  ///
  /// @param lpos is the local position to be checked
  /// @param tol1 is the tolerance applied
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const override;

  /// Return the vertices - or, the points of the extremas
  virtual const std::vector<Vector2D>
  vertices() const override;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// private helper method for inside check
  ///
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absulote tolerance parameter on the first coordinate
  /// @param tol1 absulote tolerance parameter on the second coordinate
  ///
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D& lpos, double tol0, double tol2) const;

  /// private helper method inside() method for a full symmetric trapezoid
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absulote tolerance parameter on the first coordinate
  /// @param tol1 absulote tolerance parameter on the second coordinate
  ///
  /// @return boolean indicator for the success of this operation
  bool
  insideFull(const Vector2D& lpos, double tol0 = 0., double tol1 = 0.) const;

  /// private inside() method for the triangular exclude
  /// area for an arbitrary trapezoid
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absulote tolerance parameter on the first coordinate
  /// @param tol1 absulote tolerance parameter on the second coordinate
  ///
  /// @return boolean indicator for the success of this operation
  bool
  insideExclude(const Vector2D& lpos, double tol0 = 0., double tol1 = 0.) const;

  /// private isAbove() method for checking whether a point
  /// lies above or under a straight line
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absulote tolerance parameter on the first coordinate
  /// @param tol1 absulote tolerance parameter on the second coordinate
  /// @param k is the first parameter of the parametric line equation
  /// @param d is the second parameter of the parameteric line equation
  ///
  /// @return boolean indicator for the success of this operation
  bool
  isAbove(const Vector2D& lpos, double tol0, double tol1, double k, double d)
      const;

  TDD_real_t m_alpha;  ///< private cache of angle alpha
  TDD_real_t m_beta;   ///< private cahce of angle beta
};

inline TrapezoidBounds*
TrapezoidBounds::clone() const
{
  return new TrapezoidBounds(*this);
}

inline double
TrapezoidBounds::minHalflengthX() const
{
  return m_valueStore.at(TrapezoidBounds::bv_minHalfX);
}

inline double
TrapezoidBounds::maxHalflengthX() const
{
  return m_valueStore.at(TrapezoidBounds::bv_maxHalfX);
}

inline double
TrapezoidBounds::halflengthY() const
{
  return m_valueStore.at(TrapezoidBounds::bv_halfY);
}

inline double
TrapezoidBounds::alpha() const
{
  return m_alpha;
}

inline double
TrapezoidBounds::beta() const
{
  return m_beta;
}

inline bool
TrapezoidBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  if (bcheck.bcType == 0)
    return inside(lpos, bcheck.toleranceLoc0, bcheck.toleranceLoc1);

  // a fast FALSE
  double fabsY   = fabs(lpos[Acts::eLOC_Y]);
  double max_ell = (*bcheck.lCovariance)(0, 0) > (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  double limit = bcheck.nSigmas * sqrt(max_ell);
  if (fabsY > (m_valueStore.at(TrapezoidBounds::bv_halfY) + limit))
    return false;
  // a fast FALSE
  double fabsX = fabs(lpos[Acts::eLOC_X]);
  if (fabsX > (m_valueStore.at(TrapezoidBounds::bv_maxHalfX) + limit))
    return false;
  // a fast TRUE
  double min_ell = (*bcheck.lCovariance)(0, 0) < (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  limit = bcheck.nSigmas * sqrt(min_ell);
  if (fabsX < (m_valueStore.at(TrapezoidBounds::bv_minHalfX) + limit)
      && fabsY < (m_valueStore.at(TrapezoidBounds::bv_halfY) + limit))
    return true;

  // compute KDOP and axes for surface polygon
  std::vector<KDOP>     elementKDOP(3);
  std::vector<Vector2D> elementP(4);
  float                 theta = ((*bcheck.lCovariance)(1, 0) != 0
                 && ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)) != 0)
      ? .5
          * bcheck.FastArcTan(2 * (*bcheck.lCovariance)(1, 0)
                            / ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)))
      : 0.;
  sincosCache scResult = bcheck.FastSinCos(theta);
  ActsMatrixD<2, 2> rotMatrix;
  rotMatrix << scResult.cosC, scResult.sinC, -scResult.sinC, scResult.cosC;
  ActsMatrixD<2, 2> normal;
  normal << 0, -1, 1, 0;
  // ellipse is always at (0,0), surface is moved to ellipse position and then
  // rotated
  Vector2D p;
  p << m_valueStore.at(TrapezoidBounds::bv_minHalfX),
      -m_valueStore.at(TrapezoidBounds::bv_halfY);
  elementP.at(0) = (rotMatrix * (p - lpos));
  p << -m_valueStore.at(TrapezoidBounds::bv_minHalfX),
      -m_valueStore.at(TrapezoidBounds::bv_halfY);
  elementP.at(1) = (rotMatrix * (p - lpos));
  scResult       = bcheck.FastSinCos(m_beta);
  p << m_valueStore.at(TrapezoidBounds::bv_minHalfX)
          + (2. * m_valueStore.at(TrapezoidBounds::bv_halfY))
              * (scResult.sinC / scResult.cosC),
      m_valueStore.at(TrapezoidBounds::bv_halfY);
  elementP.at(2) = (rotMatrix * (p - lpos));
  scResult       = bcheck.FastSinCos(m_alpha);
  p << -(m_valueStore.at(TrapezoidBounds::bv_minHalfX)
         + (2. * m_valueStore[TrapezoidBounds::bv_halfY])
             * (scResult.sinC / scResult.cosC)),
      m_valueStore.at(TrapezoidBounds::bv_halfY);
  elementP.at(3)             = (rotMatrix * (p - lpos));
  std::vector<Vector2D> axis = {normal * (elementP.at(1) - elementP.at(0)),
                                normal * (elementP.at(3) - elementP.at(1)),
                                normal * (elementP.at(2) - elementP.at(0))};
  bcheck.ComputeKDOP(elementP, axis, elementKDOP);
  // compute KDOP for error ellipse
  std::vector<KDOP> errelipseKDOP(3);
  bcheck.ComputeKDOP(bcheck.EllipseToPoly(3), axis, errelipseKDOP);
  // check if KDOPs overlap and return result
  return bcheck.TestKDOPKDOP(elementKDOP, errelipseKDOP);
}

inline bool
TrapezoidBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  return (fabs(lpos[Acts::eLOC_X])
          < m_valueStore.at(TrapezoidBounds::bv_maxHalfX) + tol0);
}

inline bool
TrapezoidBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  return (fabs(lpos[Acts::eLOC_Y])
          < m_valueStore.at(TrapezoidBounds::bv_halfY) + tol1);
}

inline const std::vector<Vector2D>
TrapezoidBounds::vertices() const
{
  // create the return vector
  std::vector<Vector2D> vertices;
  // fill the vertices
  vertices.reserve(4);
  vertices.push_back(
      Vector2D(m_valueStore.at(TrapezoidBounds::bv_minHalfX),
               -m_valueStore.at(TrapezoidBounds::bv_halfY)));  // [0]
  vertices.push_back(
      Vector2D(m_valueStore.at(TrapezoidBounds::bv_maxHalfX),
               m_valueStore.at(TrapezoidBounds::bv_halfY)));  // [1]
  vertices.push_back(
      Vector2D(-m_valueStore.at(TrapezoidBounds::bv_maxHalfX),
               m_valueStore.at(TrapezoidBounds::bv_halfY)));  // [1]
  vertices.push_back(
      Vector2D(-m_valueStore.at(TrapezoidBounds::bv_minHalfX),
               -m_valueStore.at(TrapezoidBounds::bv_halfY)));  // [3]
  return vertices;
}

}  // end of namespace

#endif  // ACTS_SURFACES_TRAPEZOIDBOUNDS_H
