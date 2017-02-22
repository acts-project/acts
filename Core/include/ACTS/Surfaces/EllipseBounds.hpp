// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// EllipseBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_ELLIPSEBOUNDS_H
#define ACTS_SURFACES_ELLIPSEBOUNDS_H 1

#include <cmath>
#include <cstdlib>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @class EllipseBounds
///
/// Class to describe the bounds for a planar EllipseSurface,
/// i.e. the surface between two ellipses.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phirange around the center position.
///   
/// @image html EllipseBounds.png
///  
class EllipseBounds : public PlanarBounds
{
public:
  /// @brief constants for readability
  enum BoundValues {
    bv_rMinX         = 0,
    bv_rMinY         = 1,
    bv_rMaxX         = 2,
    bv_rMaxY         = 3,
    bv_averagePhi    = 4,
    bv_halfPhiSector = 5,
    bv_length        = 6
  };

  /// Default Constructor - deleted
  EllipseBounds() = delete;

  /// Constructor for full of an ellipsoid disc
  ///
  /// @param minrad1 is the minimum radius at coorindate 1
  /// @param minrad2 is the minimum radius at coorindate 2
  /// @param maxrad1 is the minimum radius at coorindate 1
  /// @param maxrad2 is the minimum radius at coorindate 2
  /// @param avphi average phi (is set to 0. as default)
  /// @param hphisec spanning phi sector (is set to pi as default)
  EllipseBounds(double minrad1,
                double minrad2,
                double maxrad1,
                double maxrad2,
                double avphi   = 0.,
                double hphisec = M_PI);

  /// Copy constructor
  ///
  /// @param ebo is the source bounds for the copy
  EllipseBounds(const EllipseBounds& ebo)
    : PlanarBounds(ebo), m_boundingBox(0., 0.)
  {
  }

  /// Destructor
  virtual ~EllipseBounds();

  /// Assignment operator
  ///
  /// @param ebo is the source bounds for the copy
  EllipseBounds&
  operator=(const EllipseBounds& ebo);

  /// Move assignment operator
  ///
  /// @param ebo is the source bounds for the copy
  EllipseBounds&
  operator=(EllipseBounds&& ebo);

  /// Virtual constructor
  virtual EllipseBounds*
  clone() const final override;

  /// Return the type of the bounds for persistency
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Ellipse;
  }

  /// This method checks if the point given in the local coordinates is between
  /// two ellipsoids if only tol0 is given and additional in the phi sector is
  /// tol1 is given
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive

  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Check for inside first local coordinate
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol0 absolute tolerance parameter
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const final override;

  /// Check for inside second local coordinate
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol1 absolute tolerance parameter
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

  /// This method returns first inner radius
  double
  rMinX() const;

  /// This method returns second inner radius
  double
  rMinY() const;

  /// This method returns first outer radius
  double
  rMaxX() const;

  /// This method returns second outer radius
  double
  rMaxY() const;

  /// This method returns the average phi
  double
  averagePhi() const;

  /// Return the vertices - or, the points of the extremas
  virtual const std::vector<Vector2D>
  vertices() const final override;

  // Bounding box representation
  virtual const RectangleBounds&
  boundingBox() const final;

  /// This method returns the halfPhiSector which is covered by the disc
  double
  halfPhiSector() const;

  /// Output Method for std::ostream
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

private:
  /// private helper function
  ///
  /// @param lpos is the local position for checking
  /// @param tol0 is the absolute tolerance on the first parameter
  /// @param tol1 is the absolute tolerance on the second parameter
  bool
  inside(const Vector2D& lpos, double tol0, double tol1) const;

  /// helper function for squaring
  ///
  /// @param x is the input for squaring ? Why do we need this ?
  inline double
  square(double x) const
  {
    return x * x;
  };

  RectangleBounds m_boundingBox;  ///< internal bounding box cache
};

inline EllipseBounds*
EllipseBounds::clone() const
{
  return new EllipseBounds(*this);
}

inline bool
EllipseBounds::inside(const Vector2D& lpos, double tol0, double tol1) const
{
  double alpha = acos(cos(lpos[Acts::eLOC_PHI]
                          - m_valueStore.at(EllipseBounds::bv_averagePhi)));
  bool insidePhi
      = (m_valueStore.at(EllipseBounds::bv_halfPhiSector) + tol1 < M_PI)
      ? (alpha <= (m_valueStore.at(EllipseBounds::bv_halfPhiSector) + tol1))
      : 1;
  bool insideInner = (m_valueStore.at(EllipseBounds::bv_rMinX) <= tol0)
      || (m_valueStore.at(EllipseBounds::bv_rMinY) <= tol0)
      || (square(lpos[Acts::eLOC_X]
                 / (m_valueStore.at(EllipseBounds::bv_rMinX) - tol0))
              + square(lpos[Acts::eLOC_Y]
                       / (m_valueStore.at(EllipseBounds::bv_rMinY) - tol0))
          > 1);
  bool insideOuter
      = (square(lpos[Acts::eLOC_X]
                / (m_valueStore.at(EllipseBounds::bv_rMaxX) + tol0))
             + square(lpos[Acts::eLOC_Y]
                      / (m_valueStore.at(EllipseBounds::bv_rMaxY) + tol0))
         < 1);
  return (insideInner && insideOuter && insidePhi);
}

inline bool
EllipseBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  return EllipseBounds::inside(lpos, bcheck.toleranceLoc0, bcheck.toleranceLoc1);
}

inline bool
EllipseBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  bool insideInner = (m_valueStore[EllipseBounds::bv_rMinX] <= tol0)
      || (m_valueStore.at(EllipseBounds::bv_rMinY) <= tol0)
      || (square(lpos[Acts::eLOC_X]
                 / (m_valueStore.at(EllipseBounds::bv_rMinX) - tol0))
              + square(lpos[Acts::eLOC_Y]
                       / (m_valueStore.at(EllipseBounds::bv_rMinY) - tol0))
          > 1);
  bool insideOuter
      = (square(lpos[Acts::eLOC_X]
                / (m_valueStore.at(EllipseBounds::bv_rMaxX) + tol0))
             + square(lpos[Acts::eLOC_Y]
                      / (m_valueStore.at(EllipseBounds::bv_rMaxY) + tol0))
         < 1);
  return (insideInner && insideOuter);
}

inline bool
EllipseBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  double alpha = acos(cos(lpos[Acts::eLOC_PHI]
                          - m_valueStore.at(EllipseBounds::bv_averagePhi)));
  bool insidePhi
      = (m_valueStore.at(EllipseBounds::bv_halfPhiSector) + tol1 < M_PI)
      ? (alpha <= (m_valueStore.at(EllipseBounds::bv_halfPhiSector) + tol1))
      : 1;
  return insidePhi;
}

inline double
EllipseBounds::rMinX() const
{
  return m_valueStore.at(EllipseBounds::bv_rMinX);
}

inline double
EllipseBounds::rMinY() const
{
  return m_valueStore.at(EllipseBounds::bv_rMinY);
}

inline double
EllipseBounds::rMaxX() const
{
  return m_valueStore.at(EllipseBounds::bv_rMaxX);
}

inline double
EllipseBounds::rMaxY() const
{
  return m_valueStore.at(EllipseBounds::bv_rMaxY);
}

inline double
EllipseBounds::averagePhi() const
{
  return m_valueStore.at(EllipseBounds::bv_averagePhi);
}

inline double
EllipseBounds::halfPhiSector() const
{
  return m_valueStore.at(EllipseBounds::bv_halfPhiSector);
}

inline const std::vector<Vector2D>
EllipseBounds::vertices() const
{
  // create the return vector
  std::vector<Vector2D> vertices;
  // fill the vertices
  vertices.reserve(4);
  vertices.push_back(
      Vector2D(m_valueStore.at(EllipseBounds::bv_rMaxX), 0.));  // [0]
  vertices.push_back(
      Vector2D(0., m_valueStore.at(EllipseBounds::bv_rMaxY)));  // [1]
  vertices.push_back(
      Vector2D(-m_valueStore.at(EllipseBounds::bv_rMaxX), 0.));  // [2]
  vertices.push_back(
      Vector2D(0., -m_valueStore.at(EllipseBounds::bv_rMaxY)));  // [3]
  return vertices;
}

inline const RectangleBounds&
EllipseBounds::boundingBox() const
{
  return m_boundingBox;
}

}  // end of namespace

#endif  // ACTS_SURFACES_ELLIPSEBOUNDS_H
