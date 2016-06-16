// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESCONEBOUNDS_H
#define ACTS_SURFACESCONEBOUNDS_H

#include <cfloat>

#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

///  @class ConeBounds
/// 
///  Bounds for a conical Surface,
///  the opening angle is stored in \f$ \tan(\alpha) \f$ and always positively
/// defined.
///  The cone can open to both sides, steered by \f$ z_min \f$ and \f$ z_max \f$.
/// 
///  @image html ConeBounds.gif
/// 

class ConeBounds : public SurfaceBounds
{
public:
  /// @enum BoundValues for readablility
  enum BoundValues {
    bv_alpha         = 0,
    bv_minZ          = 1,
    bv_maxZ          = 2,
    bv_averagePhi    = 3,
    bv_halfPhiSector = 4,
    bv_length        = 5
  };

  /// Default Constructor is deleted
  ConeBounds() = delete;

  /// Constructor - open cone with alpha, by default a full cone
  /// but optionally can make a conical section
  /// @param alpha is the opening angle of the cone
  /// @param symm is the boolean indicating if the cone is symmetric in +/- z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param alpha is angle of the local 3D x axis (default is 0)
  ConeBounds(double alpha, bool symm, double halfphi = M_PI, double avphi = 0.);

  /// Constructor - open cone with alpha, minz and maxz, by
  /// default a full cone but can optionally make it a conical section 
  /// @param alpha is the opening angle of the cone
  /// @param zmin cone expanding from minimal z
  /// @param zmax cone expanding to maximal z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param alpha is angle of the local 3D x axis (default is 0)
  ConeBounds(double alpha,
             double zmin,
             double zmax,
             double halfphi = M_PI,
             double avphi   = 0.);

  /// Copy Constructor 
  /// @param cobo is the source bounds for the assignment             
  ConeBounds(const ConeBounds& cobo) : SurfaceBounds(cobo) {}

  /// Destructor 
  virtual ~ConeBounds();

  /// Assignment operator
  /// @param cylbo is the source bounds for the assignment             
  ConeBounds&
  operator=(const ConeBounds& cobo);

  /// Virtual constructor 
  virtual ConeBounds*
  clone() const override;

  /// Return the bounds type 
  virtual BoundsType
  type() const override
  {
    return SurfaceBounds::Cone;
  }

  /// @copydoc SurfaceBounds::inside
  virtual bool
  inside(const Vector2D&      locpo,
         const BoundaryCheck& bchk = true) const override;

  /// @copydoc SurfaceBounds::insideLoc0
  virtual bool
  insideLoc0(const Vector2D& locpo, double tol1 = 0.) const override;

  /// @copydoc SurfaceBounds::insideLoc1
  virtual bool
  insideLoc1(const Vector2D& locpo, double tol2 = 0.) const override;

  /// @copydoc SurfaceBounds::minDistance
  virtual double
  minDistance(const Vector2D& pos) const override;

  /// Return the radius at a specific z values 
  double
  r (double z) const;

  /// Return the average values for the angles (cached)
  double
  tanAlpha() const;
  double
  sinAlpha() const;
  double
  cosAlpha() const;
  double
  alpha() const;

  /// This method returns the minimum z value in the local 
  /// frame for an unbound symmetric cone, it returns -MAXBOUNDVALUE*/
  double
  minZ() const;

  /// This method returns the maximum z value in the local 
  /// frame for an unbound symmetric cone, it returns -MAXBOUNDVALUE*/
  double
  maxZ() const;

  /// This method returns the average phi value 
  /// (i.e. the "middle" phi value for the conical sector we  are describing) 
  double
  averagePhi() const;
  
  /// This method returns the half-phi width of the sector 
  /// (so that averagePhi +/- halfPhiSector gives the phi bounds of the cone)
  double
  halfPhiSector() const;

  /// Output Method for std::ostream 
  virtual std::ostream&
  dump(std::ostream& sl) const override;

private:
  std::vector<TDD_real_t> m_valueStore; ///< internal storage for the bound values
  TDD_real_t              m_tanAlpha;    ///< internal cache 
  TDD_real_t              m_sinAlpha;
  TDD_real_t              m_cosAlpha;

  /// Helper function for angle parameter initialization 
  virtual void
  initCache() override;

  /// Helpers for inside() functions 
  inline double
  minPhi() const
  {
    return m_valueStore.at(ConeBounds::bv_averagePhi)
        - m_valueStore.at(ConeBounds::bv_halfPhiSector);
  }
  inline double
  maxPhi() const
  {
    return m_valueStore.at(ConeBounds::bv_averagePhi)
        + m_valueStore.at(ConeBounds::bv_halfPhiSector);
  }
};

inline ConeBounds*
ConeBounds::clone() const
{
  return new ConeBounds(*this);
}

inline bool
ConeBounds::inside(const Vector2D& locpo, double tol1, double tol2) const
{
  double z       = locpo[Acts::eLOC_Z];
  bool   insideZ = z > (m_valueStore.at(ConeBounds::bv_minZ) - tol2)
      && z < (m_valueStore.at(ConeBounds::bv_maxZ) + tol2);
  if (!insideZ) return false;
  // TODO: Do we need some sort of "R" tolerance also here (take
  // it off the z tol2 in that case?) or does the rphi tol1 cover
  // this? (Could argue either way)
  double coneR   = z * m_tanAlpha;
  double minRPhi = coneR * minPhi() - tol1, maxRPhi = coneR * maxPhi() + tol1;
  return minRPhi < locpo[Acts::eLOC_RPHI] && locpo[Acts::eLOC_RPHI] < maxRPhi;
}

inline bool
ConeBounds::inside(const Vector3D& glopo, double tol1, double tol2) const
{
  // coords are (rphi,z)
  return inside(Vector2D(glopo.perp() * glopo.phi(), glopo.z()), tol1, tol2);
}

inline bool
ConeBounds::inside(const Vector3D& glopo, const BoundaryCheck& bchk) const
{
  // coords are (rphi,z)
  return inside(Vector2D(glopo.perp() * glopo.phi(), glopo.z()),
                bchk.toleranceLoc0,
                bchk.toleranceLoc1);
}

inline bool
ConeBounds::inside(const Vector2D& locpo, const BoundaryCheck& bchk) const
{
  return ConeBounds::inside(locpo, bchk.toleranceLoc0, bchk.toleranceLoc1);
}

inline bool
ConeBounds::insideLoc0(const Vector2D& locpo, double tol1) const
{
  double z       = locpo[Acts::eLOC_Z];
  double coneR   = z * m_tanAlpha;
  double minRPhi = coneR * minPhi() - tol1, maxRPhi = coneR * maxPhi() + tol1;
  return minRPhi < locpo[Acts::eLOC_RPHI] && locpo[Acts::eLOC_RPHI] < maxRPhi;
}

inline bool
ConeBounds::insideLoc1(const Vector2D& locpo, double tol2) const
{
  double z = locpo[Acts::eLOC_Z];
  return (z > (m_valueStore.at(ConeBounds::bv_minZ) - tol2)
          && z < (m_valueStore.at(ConeBounds::bv_maxZ) + tol2));
}

inline double
ConeBounds::r(double z) const
{
  return fabs(z * m_tanAlpha);
}

inline double
ConeBounds::tanAlpha() const
{
  return m_tanAlpha;
}

inline double
ConeBounds::sinAlpha() const
{
  return m_sinAlpha;
}

inline double
ConeBounds::cosAlpha() const
{
  return m_cosAlpha;
}

inline double
ConeBounds::alpha() const
{
  return m_valueStore.at(ConeBounds::bv_alpha);
}

inline double
ConeBounds::minZ() const
{
  return m_valueStore.at(ConeBounds::bv_minZ);
}

inline double
ConeBounds::maxZ() const
{
  return m_valueStore.at(ConeBounds::bv_maxZ);
}

inline double
ConeBounds::averagePhi() const
{
  return m_valueStore.at(ConeBounds::bv_averagePhi);
}

inline double
ConeBounds::halfPhiSector() const
{
  return m_valueStore.at(ConeBounds::bv_halfPhiSector);
}

inline void
ConeBounds::initCache()
{
  m_tanAlpha = tan(m_valueStore.at(ConeBounds::bv_alpha));
  m_sinAlpha = sin(m_valueStore.at(ConeBounds::bv_alpha));
  m_cosAlpha = cos(m_valueStore.at(ConeBounds::bv_alpha));
  // validate the halfphi
  if (m_valueStore.at(ConeBounds::bv_halfPhiSector) < 0.)
    m_valueStore.at(ConeBounds::bv_halfPhiSector)
        = -m_valueStore.at(ConeBounds::bv_halfPhiSector);
  if (m_valueStore.at(ConeBounds::bv_halfPhiSector) > M_PI)
    m_valueStore.at(ConeBounds::bv_halfPhiSector) = M_PI;
}
}

#endif  // ACTS_SURFACESCONEBOUNDS_H
