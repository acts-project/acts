// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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

  /// Default Constructor - deleted
  CylinderBounds() = delete;

  /// Constructor - full cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param halez is the half length in z
  CylinderBounds(double radius, double halez);

  /// Constructor - open cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius, double halfphi, double halez);

  /// Constructor - open cylinder
  ///
  /// @param radius is the radius of the cylinder
  /// @param avphi is the middle phi position of the segment
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius, double avphi, double halfphi, double halez);

  /// Copy Constructor
  ///
  /// @param cylbo is the source object
  CylinderBounds(const CylinderBounds& cylbo) : SurfaceBounds(cylbo) {}
  /// Destructor
  virtual ~CylinderBounds();

  /// Assignment operator
  CylinderBounds&
  operator=(const CylinderBounds& cylbo);

  /// Virtual constructor
  virtual CylinderBounds*
  clone() const final override;

  /// Return of the bounds type
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Cylinder;
  }

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Specialized method for CylinderBounds that checks if a global position
  /// is within the the cylinder cover
  ///
  /// @param pos is the position in the cylinder frame
  /// @param bcheck is the boundary check directive
  ///
  /// return boolean indicator for operation success
  bool
  inside3D(const Vector3D& pos, const BoundaryCheck& bcheck = true) const;

  /// Inside method for the second local parameter
  ///
  /// @param lpos is the local position to be checked
  /// @param tol0 is the absolute tolerance on the first parameter
  ///
  /// @return is a boolean indicating if the position is insideLoc0
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const final override;

  /// Inside method for the second local parameter
  ///
  /// @param lpos is the local position to be checked
  /// @param tol1 is the absolute tolerance on the first parameter
  ///
  /// @return is a boolean indicating if the position is insideLoc1
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const final override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const final override;

  /// This method returns the radius
  virtual double
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

  /// Output Method for std::ostream
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

private:
  /// private method for inside check
  bool
  inside(double r, double phi, double z, double tol0, double tol1) const;

  /// private method for inside check
  bool
  insideLocZ(double z, double tol1) const;

  /// private method for inside check
  bool
  inside(const Vector2D& lpos, double tol0, double tol1) const;

  /// indiciator whether to check phi or not
  bool m_checkPhi;
};

inline CylinderBounds*
CylinderBounds::clone() const
{
  return new CylinderBounds(*this);
}

inline bool
CylinderBounds::inside(const Vector2D& lpos, double tol0, double tol1) const
{
  double z       = lpos[Acts::eLOC_Z];
  bool   insideZ = insideLocZ(z, tol1);
  if (!insideZ) return false;
  // no check on Phi neccesary
  if (!m_checkPhi) return true;
  // now check insidePhi
  double localPhi
      = (lpos[Acts::eLOC_RPHI] / m_valueStore.at(CylinderBounds::bv_radius))
      - m_valueStore.at(CylinderBounds::bv_averagePhi);
  localPhi -= (localPhi > M_PI) ? 2. * M_PI : 0.;
  return (localPhi * localPhi
          < (m_valueStore.at(CylinderBounds::bv_halfPhiSector) + tol0)
              * (m_valueStore.at(CylinderBounds::bv_halfPhiSector) + tol0));
}

inline bool
CylinderBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  if (bcheck.bcType == 0 || bcheck.nSigmas == 0
      || m_valueStore.at(CylinderBounds::bv_halfPhiSector) != M_PI)
    return CylinderBounds::inside(
        lpos, bcheck.toleranceLoc0, bcheck.toleranceLoc1);

  float theta
      = ((*bcheck.lCovariance)(1, 0) != 0
         && ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)) != 0)
      ? .5 * std::atan(
                 2 * (*bcheck.lCovariance)(1, 0)
                 / ((*bcheck.lCovariance)(1, 1) - (*bcheck.lCovariance)(0, 0)))
      : 0.;
  sincosCache scResult = bcheck.FastSinCos(theta);
  double dphi    = scResult.sinC * scResult.sinC * (*bcheck.lCovariance)(0, 0);
  double dz      = scResult.cosC * scResult.cosC * (*bcheck.lCovariance)(0, 1);
  double max_ell = dphi > dz ? dphi : dz;
  double limit   = bcheck.nSigmas * sqrt(max_ell);
  return insideLocZ(lpos[Acts::eLOC_Z], limit);
}

inline bool
CylinderBounds::inside3D(const Vector3D& pos, const BoundaryCheck& bcheck) const
{
  return inside(pos.perp(),
                pos.phi(),
                pos.z(),
                bcheck.toleranceLoc0,
                bcheck.toleranceLoc0);
}

//!< @todo integrate tol0
inline bool
CylinderBounds::inside(double r,
                       double phi,
                       double z,
                       double /*tol0*/,
                       double tol1) const
{
  bool insideZ = insideLocZ(z, tol1);
  if (!insideZ) return false;
  double diffR   = (m_valueStore.at(CylinderBounds::bv_radius) - r);
  bool   insideR = diffR * diffR < s_onSurfaceTolerance * s_onSurfaceTolerance;
  if (!insideR) return false;
  // now check insidePhi if needed
  if (!m_checkPhi) return true;
  // phi needs to be checked
  double localPhi = phi - m_valueStore.at(CylinderBounds::bv_averagePhi);
  localPhi -= (localPhi > M_PI) ? 2. * M_PI : 0.;
  return (localPhi * localPhi
          < m_valueStore.at(CylinderBounds::bv_halfPhiSector)
              * m_valueStore.at(CylinderBounds::bv_halfPhiSector));
}

inline bool
CylinderBounds::insideLocZ(double z, double tol1) const
{
  return (m_valueStore.at(CylinderBounds::bv_halfZ) + tol1) - std::abs(z) > 0.;
}

inline bool
CylinderBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  bool insideRphi = false;
  if (std::abs(m_valueStore.at(CylinderBounds::bv_averagePhi)) < 10e-7)
    insideRphi = (std::abs(lpos[Acts::eLOC_RPHI]
                           / m_valueStore.at(CylinderBounds::bv_radius))
                  < (m_valueStore.at(CylinderBounds::bv_halfPhiSector) + tol0));
  else {
    double localPhi
        = (lpos[Acts::eLOC_RPHI] / m_valueStore.at(CylinderBounds::bv_radius))
        - m_valueStore.at(CylinderBounds::bv_averagePhi);
    localPhi -= (localPhi > M_PI) ? 2. * M_PI : 0.;
    insideRphi = (localPhi
                  < (m_valueStore.at(CylinderBounds::bv_halfPhiSector) + tol0));
  }
  return (insideRphi);
}

inline bool
CylinderBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  return insideLocZ(lpos[Acts::eLOC_Z], tol1);
}

inline double
CylinderBounds::r() const
{
  return m_valueStore.at(CylinderBounds::bv_radius);
}

inline double
CylinderBounds::averagePhi() const
{
  return m_valueStore.at(CylinderBounds::bv_averagePhi);
}

inline double
CylinderBounds::halfPhiSector() const
{
  return m_valueStore.at(CylinderBounds::bv_halfPhiSector);
}

inline double
CylinderBounds::halflengthZ() const
{
  return m_valueStore.at(CylinderBounds::bv_halfZ);
}
}

#endif  // ACTS_SURFACES_CYLINDERBOUNDS_H
