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
  /// @param radius is the radius of the cylinder
  /// @param halez is the half length in z 
  CylinderBounds(double radius, double halez);

  /// Constructor - open cylinder
  /// @param radius is the radius of the cylinder
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius, double halfphi, double halez);

  /// Constructor - open cylinder
  /// @param radius is the radius of the cylinder
  /// @param avphi is the middle phi position of the segment
  /// @param halfphi is the half opening angle
  /// @param halez is the half length in z
  CylinderBounds(double radius, double avphi, double halfphi, double halez);

  /// Copy Constructor 
  CylinderBounds(const CylinderBounds& cylbo) : SurfaceBounds(cylbo) {}

  /// Destructor
  virtual ~CylinderBounds();

  /// Assignment operator
  CylinderBounds&
  operator=(const CylinderBounds& cylbo);

  /// Virtual constructor
  virtual CylinderBounds*
  clone() const override;

  /// Return of the bounds type
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Cylinder;
  }

  /// @copydoc SurfaceBounds::inside
  inside(const Vector2D& lpos, const BoundaryCheck& bchk) const override;

  /// specialized method for CylinderBounds 
  bool
  inside3D(const Vector3D& gp, const BoundaryCheck& bchk = true) const;

  /// @copydoc Surface::insideLoc0
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const override;

  /// @copydoc Surface::insideLoc1
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const override;

  /// Minimal distance to boundary 
  /// return minimal distance to boundary ( > 0 if outside and <=0 if inside) 
  virtual double
  minDistance(const Vector2D& pos) const override;
  
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
  dump(std::ostream& sl) const override;

private:
  /// private method for inside check
  bool
  inside(double r, double phi, double z, double tol0, double tol1) const;
  
  /// private method for inside check
  bool
  insideLocZ(double z, double tol1) const;
  
  /// private method for inside check
  bool inside(const Vector2D& lpos, double tol0, double tol1) const;
  
  /// indiciator whether to check phi or not
  bool                    m_checkPhi;
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
CylinderBounds::inside(const Vector2D& lpos, const BoundaryCheck& bchk) const
{
  if (bchk.bcType == 0 || bchk.nSigmas == 0
      || m_valueStore.at(CylinderBounds::bv_halfPhiSector) != M_PI)
    return CylinderBounds::inside(
        lpos, bchk.toleranceLoc0, bchk.toleranceLoc1);

  float theta = (bchk.lCovariance(1, 0) != 0
                 && (bchk.lCovariance(1, 1) - bchk.lCovariance(0, 0)) != 0)
      ? .5
          * bchk.FastArcTan(2 * bchk.lCovariance(1, 0)
                            / (bchk.lCovariance(1, 1) - bchk.lCovariance(0, 0)))
      : 0.;
  sincosCache scResult = bchk.FastSinCos(theta);
  double      dphi     = scResult.sinC * scResult.sinC * bchk.lCovariance(0, 0);
  double      dz       = scResult.cosC * scResult.cosC * bchk.lCovariance(0, 1);
  double      max_ell  = dphi > dz ? dphi : dz;
  double      limit    = bchk.nSigmas * sqrt(max_ell);
  return insideLocZ(lpos[Acts::eLOC_Z], limit);
}

inline bool
CylinderBounds::inside3D(const Vector3D& glopo, const BoundaryCheck& bchk) const
{
  return inside(glopo.perp(),
                glopo.phi(),
                glopo.z(),
                bchk.toleranceLoc0,
                bchk.toleranceLoc0);
}

//!< @TODO integrate tol0
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
  return (m_valueStore.at(CylinderBounds::bv_halfZ) + tol1) - fabs(z) > 0.;
}

inline bool
CylinderBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  bool insideRphi = false;
  if (fabs(m_valueStore.at(CylinderBounds::bv_averagePhi)) < 10e-7)
    insideRphi
        = (fabs(lpos[Acts::eLOC_RPHI]
                / m_valueStore.at(CylinderBounds::bv_radius))
           < (m_valueStore.at(CylinderBounds::bv_halfPhiSector) + tol0));
  else {
    double localPhi
        = (lpos[Acts::eLOC_RPHI] / m_valueStore.at(CylinderBounds::bv_radius))
        - m_valueStore.at(CylinderBounds::bv_averagePhi);
    localPhi -= (localPhi > M_PI) ? 2. * M_PI : 0.;
    insideRphi = (localPhi < (m_valueStore.at(CylinderBounds::bv_halfPhiSector)
                              + tol0));
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
