// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_SUBTRACTEDVOLUMEBOUNDS_H
#define ACTS_VOLUMES_SUBTRACTEDVOLUMEBOUNDS_H 1

// Geometry module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/Volume.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
// Core module

namespace Acts {

class SurfaceBounds;
class Volume;
class Surface;

/**
 @class SubtractedVolumeBounds

 Bounds for a generic subtracted volume, the decomposeToSurfaces method creates
 a
 vector of n surfaces (n1+n2-n_subtracted):

  BoundarySurfaceFace [index]: [n1] surfaces from 'outer' volume
                               [n1+n2-n_subtr] remaining surfaces (after
 subtraction) from 'inner' volume

  */

class SubtractedVolumeBounds : public VolumeBounds
{
public:
  /**Default Constructor*/
  SubtractedVolumeBounds();

  /**Constructor - the box boundaries */
  SubtractedVolumeBounds(Volume* outerVol, Volume* innerVol);

  /**Copy Constructor */
  SubtractedVolumeBounds(const SubtractedVolumeBounds& bobo);

  /**Destructor */
  virtual ~SubtractedVolumeBounds();

  /**Assignment operator*/
  SubtractedVolumeBounds&
  operator=(const SubtractedVolumeBounds& bobo);

  /**Virtual constructor */
  SubtractedVolumeBounds*
  clone() const override;

  /**This method checks if position in the 3D volume frame is inside the
   * cylinder*/
  bool
  inside(const Vector3D&, double tol = 0.) const override;

  /** Method to decompose the Bounds into boundarySurfaces */
  const std::vector<const Acts::Surface*>*
  decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

  /**This method returns the outer Volume*/
  Volume*
  outer() const;

  /**This method returns the inner Volume*/
  Volume*
  inner() const;

  /**This method returns bounds orientation*/
  const std::vector<bool>
  boundsOrientation() const;

  /** Output Method for std::ostream */
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  Acts::Volume*
  createSubtractedVolume(const Transform3D& transf,
                         Acts::Volume*      subtrVol) const;

  Volume* m_outer;
  Volume* m_inner;

  mutable std::vector<bool> m_boundsOrientation;
};

inline SubtractedVolumeBounds*
SubtractedVolumeBounds::clone() const
{
  return new SubtractedVolumeBounds(*this);
}

inline bool
SubtractedVolumeBounds::inside(const Vector3D& pos, double tol) const
{
  return (m_outer->inside(pos, tol) && !m_inner->inside(pos, -tol));
}

inline Volume*
SubtractedVolumeBounds::outer() const
{
  return m_outer;
}

inline Volume*
SubtractedVolumeBounds::inner() const
{
  return m_inner;
}

inline const std::vector<bool>
SubtractedVolumeBounds::boundsOrientation() const
{
  return (m_boundsOrientation);
}
}

#endif  // ACTS_VOLUMES_SUBTRACTEDVOLUMEBOUNDS_H
