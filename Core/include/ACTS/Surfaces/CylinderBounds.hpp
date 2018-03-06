// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_CYLINDERBOUNDS_H
#define ACTS_SURFACES_CYLINDERBOUNDS_H 1

#include <cmath>

#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/detail/periodic.hpp"

namespace Acts {

/// @class CylinderBounds
///
/// Bounds for a cylindrical Surface.
///
/// These bounds may be used for a CylinderSurface
/// In case of bounds for a StraightLineSurface the radius determines the radius
/// within a localPosition
/// is regarded as inside bounds.
///
/// CylinderBounds also enhance the possibility of a cylinder segment with an
/// opening angle @f$ 2\cdot\phi_{half}@f$
/// around an average @f$ \phi @f$ angle @f$ \phi_{ave} @f$.
///
/// @todo update the documentation picture for cylinder segments
///
/// @image html CylinderBounds.gif

class CylinderBounds : public SurfaceBounds
{
public:
  /// @enum BoundValues for readablility
  /// nested enumeration object
  enum BoundValues {
    bv_radius        = 0,
    bv_averagePhi    = 1,
    bv_halfPhiSector = 2,
    bv_halfZ         = 3,
    bv_length        = 4
  };

  CylinderBounds() = delete;

  /// Constructor - full cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param halez is the half length in z
  CylinderBounds(double radius, double halfZ);

  /// Constructor - open cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius, double halfPhi, double haleZ);

  /// Constructor - open cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param avphi is the middle phi position of the segment
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius,
                 double averagePhi,
                 double halfPhi,
                 double halfZ);

  CylinderBounds(const variant_data& data);

  virtual ~CylinderBounds();

  virtual CylinderBounds*
  clone() const final override;

  virtual BoundsType
  type() const final override;

  virtual std::vector<TDD_real_t>
  valueStore() const final override;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Specialized method for CylinderBounds that checks if a global position
  /// is within the the cylinder cover
  ///
  /// @param pos is the position in the cylinder frame
  /// @param bcheck is the boundary check directive
  /// @return boolean indicator for operation success
  bool
  inside3D(const Vector3D& pos, const BoundaryCheck& bcheck = true) const;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const final override;

  /// Output Method for std::ostream
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

  /// This method returns the radius
  double
  r() const;

  /// This method returns the average phi
  double
  averagePhi() const;

  /// This method returns the halfPhiSector angle
  double
  halfPhiSector() const;

  /// This method returns the halflengthZ
  double
  halflengthZ() const;
  
  variant_data
  toVariantData() const;

private:
  double m_radius, m_avgPhi, m_halfPhi, m_halfZ;

  Vector2D
  shifted(const Vector2D& lpos) const;
  ActsSymMatrixD<2>
  jacobian() const;
};

inline double
CylinderBounds::r() const
{
  return m_radius;
}

inline double
CylinderBounds::averagePhi() const
{
  return m_avgPhi;
}

inline double
CylinderBounds::halfPhiSector() const
{
  return m_halfPhi;
}

inline double
CylinderBounds::halflengthZ() const
{
  return m_halfZ;
}
}

#endif  // ACTS_SURFACES_CYLINDERBOUNDS_H
