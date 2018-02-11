// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
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
#include "ACTS/Utilities/VariantDataFwd.hpp"

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

  EllipseBounds() = delete;

  /// Constructor for full of an ellipsoid disc
  ///
  /// @param minRadius0 is the minimum radius along coordinate 0
  /// @param minRadius1 is the minimum radius along coorindate 1
  /// @param maxRadius0 is the minimum radius at coorindate 0
  /// @param maxRadius1 is the minimum radius at coorindate 1
  /// @param averagePhi average phi (is set to 0. as default)
  /// @param halfPhi    spanning phi sector (is set to pi as default)
  EllipseBounds(double minRadius0,
                double minRadius1,
                double maxRadius0,
                double maxRadius1,
                double averagePhi = 0.,
                double halfPhi    = M_PI);

  virtual ~EllipseBounds();

  virtual EllipseBounds*
  clone() const final override;

  virtual BoundsType
  type() const final override;

  virtual std::vector<TDD_real_t>
  valueStore() const final override;

  /// This method checks if the point given in the local coordinates is between
  /// two ellipsoids if only tol0 is given and additional in the phi sector is
  /// tol1 is given
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const final override;

  /// Return the vertices - or, the points of the extremas
  virtual std::vector<Vector2D>
  vertices() const final override;

  // Bounding box representation
  virtual const RectangleBounds&
  boundingBox() const final;

  /// Output Method for std::ostream
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

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

  /// This method returns the halfPhiSector which is covered by the disc
  double
  halfPhiSector() const;
  
  variant_data
  toVariantData() const override;

private:
  double          m_rMinX, m_rMinY, m_rMaxX, m_rMaxY, m_avgPhi, m_halfPhi;
  RectangleBounds m_boundingBox;
};

inline double
EllipseBounds::rMinX() const
{
  return m_rMinX;
}

inline double
EllipseBounds::rMinY() const
{
  return m_rMinY;
}

inline double
EllipseBounds::rMaxX() const
{
  return m_rMaxX;
}

inline double
EllipseBounds::rMaxY() const
{
  return m_rMaxY;
}

inline double
EllipseBounds::averagePhi() const
{
  return m_avgPhi;
}

inline double
EllipseBounds::halfPhiSector() const
{
  return m_halfPhi;
}

}  // end of namespace

#endif  // ACTS_SURFACES_ELLIPSEBOUNDS_H
