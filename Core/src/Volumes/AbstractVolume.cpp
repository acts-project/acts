// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AbstractVolume.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Volumes/BoundaryCylinderSurface.hpp"
#include "ACTS/Volumes/BoundaryDiscSurface.hpp"
#include "ACTS/Volumes/BoundaryPlaneSurface.hpp"
#include "ACTS/Volumes/BoundarySurface.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
#include <iostream>

Acts::AbstractVolume::AbstractVolume(
    std::shared_ptr<Transform3D>        htrans,
    std::shared_ptr<const VolumeBounds> volbounds)
  : Volume(htrans, volbounds)
{
  createBoundarySurfaces();
}

Acts::AbstractVolume::~AbstractVolume()
{
}


const std::vector< BoundarySurfacePtr >&
Acts::AbstractVolume::boundarySurfaces() const
{
  return m_boundarySurfaces;
}

void
Acts::AbstractVolume::createBoundarySurfaces()
{
  // transform Surfaces To BoundarySurfaces
  const std::vector<const Surface*>* surfaces
      = Volume::volumeBounds().decomposeToSurfaces(m_transform);
  
  std::vector<const Surface*>::const_iterator surfIter
      = surfaces->begin();

  // counter to flip the inner/outer position for Cylinders
  int sfCounter = 0;
  int sfNumber  = surfaces->size();

  for (; surfIter != surfaces->end(); ++surfIter) {
    sfCounter++;
    const PlaneSurface* psf
        = dynamic_cast<const PlaneSurface*>(*surfIter);
    if (psf) {
      m_boundarySurfaces->push_back(
          BoundarySurfacePtr(
              new BoundaryPlaneSurface<AbstractVolume>(
                  this, 0, *psf)));
      delete psf;
      continue;
    }
    const DiscSurface* dsf
        = dynamic_cast<const DiscSurface*>(*surfIter);
    if (dsf) {
      m_boundarySurfaces->push_back(
          BoundarySurfacePtr(
              new BoundaryDiscSurface<AbstractVolume>(
                  this, 0, *dsf)));
      delete dsf;
      continue;
    }
    const CylinderSurface* csf
        = dynamic_cast<const CylinderSurface*>(*surfIter);
    if (csf) {
      AbstractVolume* inner = (sfCounter == 3 && sfNumber > 3) ? 0 : this;
      AbstractVolume* outer = (inner) ? 0 : this;
      m_boundarySurfaces->push_back(
          BoundarySurfacePtr(
              new BoundaryCylinderSurface<AbstractVolume>(
                  inner, outer, *csf)));
      delete csf;
      continue;
    }
  }

  delete surfaces;
}
