// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SlidingCylinderSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SLIDINGCYLINDERSURFACE_H
#define ACTS_SURFACES_SLIDINGCYLINDERSURFACE_H 1

#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

class Identifier;

namespace Acts {

/**
 @class SlidingCylinderSurface

 Class for a Calo CylinderSurface with variable depth in the ATLAS detector.
 The variable depth is stored as a binned vector of radial corrections.
 Local eta bin is defined by base curvature and z position in base transform (
 corrected for misalignement ).
 It inherits from CylinderSurface.

 */

class SlidingCylinderSurface : public CylinderSurface
{
public:
  /** Default Constructor - needed for persistency*/
  SlidingCylinderSurface();

  /** Copy Constructor*/
  SlidingCylinderSurface(const SlidingCylinderSurface& psf);

  /** Copy Constructor with shift*/
  SlidingCylinderSurface(const SlidingCylinderSurface& psf,
                         const Transform3D&            transf);

  /** Constructor */
  SlidingCylinderSurface(const CylinderSurface&    surf,
                         BinUtility*               bu     = nullptr,
                         const std::vector<float>* offset = nullptr,
                         Transform3D*              align  = nullptr);

  /** Destructor*/
  virtual ~SlidingCylinderSurface();

  /** Virtual constructor - shift can be optionally applied */
  virtual SlidingCylinderSurface*
  clone(const Transform3D* shift = nullptr) const;

  /** Assignment operator*/
  SlidingCylinderSurface&
  operator=(const SlidingCylinderSurface& psf);

  /** Equality operator*/
  bool
  operator==(const Surface& sf) const;

  /** This method returns true if the GlobalPosition is on the Surface for both,
    within
    or without check of whether the local position is inside boundaries or not
    */
  bool
  isOnSurface(const Vector3D& glopo, const BoundaryCheck& bchk = true) const;

  /** Specialized for DiscSurface: LocalToGlobal method without dynamic memory
   * allocation */
  void
  localToGlobal(const Vector2D& locp,
                const Vector3D& mom,
                Vector3D&       glob) const;

  /** Specialized for DiscSurface: GlobalToLocal method without dynamic memory
   * allocation - boolean checks if on surface */
  bool
  globalToLocal(const Vector3D& glob, const Vector3D& mom, Vector2D& loc) const;

  /**This method allows access to the bin utility*/
  const BinUtility*
  binUtility() const
  {
    return m_etaBin;
  }

  /**This method allows access to the radial offset values*/
  const std::vector<float>*
  offset() const
  {
    return m_depth;
  }

  /** Return properly formatted class name for screen output */
  std::string
  name() const
  {
    return "Acts::SlidingCylinderSurface";
  }

protected:
  const std::vector<float>* m_depth;
  BinUtility*               m_etaBin;
  Transform3D*              m_align;
};

inline SlidingCylinderSurface*
SlidingCylinderSurface::clone(const Transform3D* shift) const
{
  if (shift) new SlidingCylinderSurface(*this, *shift);
  return new SlidingCylinderSurface(*this);
}

}  // end of namespace

#endif  // ACTS_SURFACES_SLIDINGCYLINDERSURFACE_H
