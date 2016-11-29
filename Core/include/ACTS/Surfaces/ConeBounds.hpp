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
///  The cone can open to both sides, steered by \f$ z_min \f$ and \f$ z_max
///  \f$.
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
  ///
  /// @param alpha is the opening angle of the cone
  /// @param symm is the boolean indicating if the cone is symmetric in +/- z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, bool symm, double halfphi = M_PI, double avphi = 0.);

  /// Constructor - open cone with alpha, minz and maxz, by
  /// default a full cone but can optionally make it a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param zmin cone expanding from minimal z
  /// @param zmax cone expanding to maximal z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha,
             double zmin,
             double zmax,
             double halfphi = M_PI,
             double avphi   = 0.);

  /// Copy Constructor
  ///
  /// @param cobo is the source bounds for the assignment
  ConeBounds(const ConeBounds& cobo)
    : SurfaceBounds(cobo)
    , m_tanAlpha(cobo.m_tanAlpha)
    , m_sinAlpha(cobo.m_sinAlpha)
    , m_cosAlpha(cobo.m_cosAlpha)
  {
  }

  /// Destructor
  virtual ~ConeBounds();

  /// Assignment operator
  /// @param cobo is the source bounds for the assignment
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

  /// inside method for local position
  ///
  /// @param lpos is the local position to be checked
  /// @param bcheck is the boundary check directive
  ///
  /// @return is a boolean indicating if the position is inside
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck = true) const override;

  /// Inside method for the first local parameter
  ///
  /// @param lpos is the local position to be checked
  /// @param tol0 is the absolute tolerance on the first parameter
  ///
  /// @return is a boolean indicating if the position is insideLoc0
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const override;

  /// Inside method for the second local parameter
  ///
  /// @param lpos is the local position to be checked
  /// @param tol1 is the absolute tolerance on the first parameter
  ///
  /// @return is a boolean indicating if the position is insideLoc1
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const override;

  /// Return the radius at a specific z values
  ///
  /// @param z is the z value for which r is requested
  ///
  /// @return is the r value associated with z
  double
  r(double z) const;

  /// Return the average values for the angles (cached)
  double
  tanAlpha() const;

  /// Return the average values for the angles (cached)
  double
  sinAlpha() const;

  /// Return the average values for the angles (cached)
  double
  cosAlpha() const;

  /// Return the average values for the angles (cached)
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
  ///
  /// @param sl is the ostrea into which the dump is done
  ///
  /// @return is the input obect
  virtual std::ostream&
  dump(std::ostream& sl) const override;

private:
  /// private helper method
  bool
  inside(const Vector2D& lpos, double tol0, double tol1) const;

  std::vector<TDD_real_t>
             m_valueStore;  ///< internal storage for the bound values
  TDD_real_t m_tanAlpha;    ///< internal cache
  TDD_real_t m_sinAlpha;    ///< internal cache
  TDD_real_t m_cosAlpha;    ///< internal cache

  /// Helper function for angle parameter initialization
  virtual void
  initCache();

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
ConeBounds::inside(const Vector2D& lpos, double tol0, double tol1) const
{
  double z       = lpos[Acts::eLOC_Z];
  bool   insideZ = z > (m_valueStore.at(ConeBounds::bv_minZ) - tol1)
      && z < (m_valueStore.at(ConeBounds::bv_maxZ) + tol1);
  if (!insideZ) return false;
  // TODO: Do we need some sort of "R" tolerance also here (take
  // it off the z tol1 in that case?) or does the rphi tol0 cover
  // this? (Could argue either way)
  double coneR   = z * m_tanAlpha;
  double minRPhi = coneR * minPhi() - tol0, maxRPhi = coneR * maxPhi() + tol0;
  return minRPhi < lpos[Acts::eLOC_RPHI] && lpos[Acts::eLOC_RPHI] < maxRPhi;
}

inline bool
ConeBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  return ConeBounds::inside(lpos, bcheck.toleranceLoc0, bcheck.toleranceLoc1);
}

inline bool
ConeBounds::insideLoc0(const Vector2D& lpos, double tol0) const
{
  double z       = lpos[Acts::eLOC_Z];
  double coneR   = z * m_tanAlpha;
  double minRPhi = coneR * minPhi() - tol0, maxRPhi = coneR * maxPhi() + tol0;
  return minRPhi < lpos[Acts::eLOC_RPHI] && lpos[Acts::eLOC_RPHI] < maxRPhi;
}

inline bool
ConeBounds::insideLoc1(const Vector2D& lpos, double tol1) const
{
  double z = lpos[Acts::eLOC_Z];
  return (z > (m_valueStore.at(ConeBounds::bv_minZ) - tol1)
          && z < (m_valueStore.at(ConeBounds::bv_maxZ) + tol1));
}

inline double
ConeBounds::r(double z) const
{
  return std::abs(z * m_tanAlpha);
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
