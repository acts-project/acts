// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiamondBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_DIAMONDDBOUNDS_H
#define ACTS_SURFACES_DIAMONDDBOUNDS_H 1

#include <math.h>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

///
/// @class DiamondBounds
///  
/// Bounds for a double trapezoidal ("diamond"), planar Surface.
///
class DiamondBounds : public PlanarBounds
{
public:
  /// @enum BoundValues for better readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_medHalfX = 1,
    bv_maxHalfX = 2,
    bv_halfY1   = 3,
    bv_halfY2   = 4,
    bv_length   = 5
  };

  /// Default Constructor is deleted
  DiamondBounds() = delete;

  /// Constructor for symmetric Diamond
  ///
  /// @param minhalex is the halflength in x at minimal y
  /// @param medhalex is the halflength in x at y = 0
  /// @param maxhalex is the halflength in x at maximal y
  /// @param haley1 is the halflength into y < 0
  /// @param haley2 is the halflength into y > 0
  DiamondBounds(double minhalex,
                double medhalex,
                double maxhalex,
                double haley1,
                double haley2);

  /// Copy constructor
  ///
  /// @param diabo are the source bounds for the copy
  DiamondBounds(const DiamondBounds& diabo);

  /// Destructor
  virtual ~DiamondBounds();

  /// Virtual constructor
  DiamondBounds*
  clone() const override;

  /// Assignment operator
  ///
  /// @param diabo are the source bounds for the copy
  DiamondBounds&
  operator=(const DiamondBounds& diabo);

  /// Comparison (Equality) operator
  ///
  /// @param sbo are the source bounds for check
  virtual bool
  operator==(const SurfaceBounds& sbo) const override;

  /// Return the bounds type
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Diamond;
  }

  /// This method returns the halflength in X at minimal Y
  /// (first coordinate of local surface frame)
  double
  minHalflengthX() const;

  /// This method returns the (maximal) halflength in X
  /// (first coordinate of local surface frame)
  double
  medHalflengthX() const;

  /// This method returns the halflength in X at maximal Y
  /// (first coordinate of local surface frame)
  double
  maxHalflengthX() const;

  /// This method returns the halflength in Y of trapezoid at negative Y
  double
  halflengthY1() const;

  /// This method returns the halflength in Y of trapezoid at positive Y
  double
  halflengthY2() const;

  /// This method returns the opening angle alpha in point A
  double
  alpha1() const;

  /// This method returns the opening angle alpha in point A'
  double
  alpha2() const;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const override;

  ///  This method checks inside bounds in loc0
  /// - loc0/loc1 correspond to the natural coordinates of the surface
  /// - As loc0/loc1 are correlated the single check doesn't make sense :
  /// -> check is done on enclosing Rectangle !
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 is the absolute tolerance
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const override;

  ///  This method checks inside bounds in loc1
  /// - loc0/loc1 correspond to the natural coordinates of the surface
  /// - As loc0/loc1 are correlated the single check doesn't make sense :
  /// -> check is done on enclosing Rectangle !
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol1 is the absolute tolerance
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
  /// @param sl is the ostream in which it is dumped
  virtual std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// private helper method
  bool
  inside(const Vector2D& lpos, double tol0 = 0., double tol1 = 0.) const;

  /// inside() method for a full symmetric diamond
  bool
  insideFull(const Vector2D& lpos, double tol0 = 0., double tol1 = 0.) const;

  /// initialize the alpha1/2 cache - needed also for object persistency
  virtual void
  initCache();

  std::vector<TDD_real_t> m_valueStore;  ///< internal parameter store
  TDD_real_t              m_alpha1;  ///< internal parameter cache for alpha1
  TDD_real_t              m_alpha2;  ///< internal parameter cahce for alpha2
};

inline DiamondBounds*
DiamondBounds::clone() const
{
  return new DiamondBounds(*this);
}

inline double
DiamondBounds::minHalflengthX() const
{
  return m_valueStore.at(DiamondBounds::bv_minHalfX);
}

inline double
DiamondBounds::medHalflengthX() const
{
  return m_valueStore.at(DiamondBounds::bv_medHalfX);
}

inline double
DiamondBounds::maxHalflengthX() const
{
  return m_valueStore.at(DiamondBounds::bv_maxHalfX);
}

inline double
DiamondBounds::halflengthY1() const
{
  return m_valueStore.at(DiamondBounds::bv_halfY1);
}

inline double
DiamondBounds::halflengthY2() const
{
  return m_valueStore.at(DiamondBounds::bv_halfY2);
}

inline bool
DiamondBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  if (bcheck.bcType == 0)
    return DiamondBounds::inside(lpos, bcheck.toleranceLoc0, bcheck.toleranceLoc1);

  // a fast FALSE
  double max_ell = (*bcheck.lCovariance)(0, 0) > (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  double limit = bcheck.nSigmas * sqrt(max_ell);
  if (lpos[Acts::eLOC_Y]
      < -2 * m_valueStore.at(DiamondBounds::bv_halfY1) - limit)
    return false;
  if (lpos[Acts::eLOC_Y]
      > 2 * m_valueStore.at(DiamondBounds::bv_halfY2) + limit)
    return false;
  // a fast FALSE
  double fabsX = fabs(lpos[Acts::eLOC_X]);
  if (fabsX > (m_valueStore.at(DiamondBounds::bv_medHalfX) + limit))
    return false;
  // a fast TRUE
  double min_ell = (*bcheck.lCovariance)(0, 0) < (*bcheck.lCovariance)(1, 1)
      ? (*bcheck.lCovariance)(0, 0)
      : (*bcheck.lCovariance)(1, 1);
  limit = bcheck.nSigmas * sqrt(min_ell);
  if (fabsX < (fmin(m_valueStore.at(DiamondBounds::bv_minHalfX),
                    m_valueStore.at(DiamondBounds::bv_maxHalfX))
               - limit))
    return true;
  // a fast TRUE
  if (fabs(lpos[Acts::eLOC_Y])
      < (fmin(m_valueStore.at(DiamondBounds::bv_halfY1),
              m_valueStore.at(DiamondBounds::bv_halfY2))
         - limit))
    return true;

  // compute KDOP and axes for surface polygon
  std::vector<KDOP>     elementKDOP(5);
  std::vector<Vector2D> elementP(6);
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
  p << -m_valueStore.at(DiamondBounds::bv_minHalfX),
      -2. * m_valueStore.at(DiamondBounds::bv_halfY1);
  elementP.at(0) = (rotMatrix * (p - lpos));
  p << -m_valueStore.at(DiamondBounds::bv_medHalfX), 0.;
  elementP.at(1) = (rotMatrix * (p - lpos));
  p << -m_valueStore.at(DiamondBounds::bv_maxHalfX),
      2. * m_valueStore.at(DiamondBounds::bv_halfY2);
  elementP.at(2) = (rotMatrix * (p - lpos));
  p << m_valueStore.at(DiamondBounds::bv_maxHalfX),
      2. * m_valueStore.at(DiamondBounds::bv_halfY2);
  elementP.at(3) = (rotMatrix * (p - lpos));
  p << m_valueStore.at(DiamondBounds::bv_medHalfX), 0.;
  elementP.at(4) = (rotMatrix * (p - lpos));
  p << m_valueStore.at(DiamondBounds::bv_minHalfX),
      -2. * m_valueStore.at(DiamondBounds::bv_halfY1);
  elementP.at(5)             = (rotMatrix * (p - lpos));
  std::vector<Vector2D> axis = {normal * (elementP.at(1) - elementP.at(0)),
                                normal * (elementP.at(2) - elementP.at(1)),
                                normal * (elementP.at(3) - elementP.at(2)),
                                normal * (elementP.at(4) - elementP.at(3)),
                                normal * (elementP.at(5) - elementP.at(4))};
  bcheck.ComputeKDOP(elementP, axis, elementKDOP);
  // compute KDOP for error ellipse
  std::vector<KDOP> errelipseKDOP(5);
  bcheck.ComputeKDOP(bcheck.EllipseToPoly(3), axis, errelipseKDOP);
  // check if KDOPs overlap and return result
  return bcheck.TestKDOPKDOP(elementKDOP, errelipseKDOP);
}

inline bool
DiamondBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  return (fabs(lpos[Acts::eLOC_X])
          < m_valueStore.at(DiamondBounds::bv_medHalfX) + tol0);
}

inline bool
DiamondBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  return ((lpos[Acts::eLOC_Y]
           > -2. * m_valueStore.at(DiamondBounds::bv_halfY1) - tol1)
          && (lpos[Acts::eLOC_Y]
              < 2. * m_valueStore.at(DiamondBounds::bv_halfY2) + tol1));
}

inline void
DiamondBounds::initCache()
{
  m_alpha1 = atan2(m_valueStore.at(DiamondBounds::bv_medHalfX)
                       - m_valueStore.at(DiamondBounds::bv_minHalfX),
                   2. * m_valueStore.at(DiamondBounds::bv_halfY1));
  m_alpha2 = atan2(m_valueStore.at(DiamondBounds::bv_medHalfX)
                       - m_valueStore.at(DiamondBounds::bv_maxHalfX),
                   2. * m_valueStore.at(DiamondBounds::bv_halfY2));
}

inline const std::vector<Vector2D>
DiamondBounds::vertices() const
{
  // create the return vector
  std::vector<Vector2D> vertices;
  // fill the vertices
  vertices.reserve(6);
  vertices.push_back(
      Vector2D(m_valueStore.at(DiamondBounds::bv_minHalfX),
               -m_valueStore.at(DiamondBounds::bv_halfY1)));  // [0]
  vertices.push_back(
      Vector2D(m_valueStore.at(DiamondBounds::bv_medHalfX), 0.));  // [1]
  vertices.push_back(
      Vector2D(m_valueStore.at(DiamondBounds::bv_maxHalfX),
               m_valueStore.at(DiamondBounds::bv_halfY2)));  // [2]
  vertices.push_back(
      Vector2D(-m_valueStore.at(DiamondBounds::bv_maxHalfX),
               m_valueStore.at(DiamondBounds::bv_halfY2)));  // [3]
  vertices.push_back(
      Vector2D(-m_valueStore.at(DiamondBounds::bv_medHalfX), 0.));  // [4]
  vertices.push_back(
      Vector2D(-m_valueStore.at(DiamondBounds::bv_minHalfX),
               -m_valueStore.at(DiamondBounds::bv_halfY1)));  // [5]
  return vertices;
}

}  // end of namespace

#endif  // ACTS_SURFACES_DIAMONDBOUNDS_H
