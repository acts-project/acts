// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AbstractVolume.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "Acts/Volumes/AbstractVolume.hpp"
#include <iostream>
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

Acts::AbstractVolume::AbstractVolume(
    std::shared_ptr<const Transform3D>  htrans,
    std::shared_ptr<const VolumeBounds> volbounds)
  : Volume(htrans, volbounds)
{
  createBoundarySurfaces();
}

Acts::AbstractVolume::~AbstractVolume()
{
}

const std::vector<Acts::BoundarySurfacePtr>&
Acts::AbstractVolume::boundarySurfaces() const
{
  return m_boundarySurfaces;
}

void
Acts::AbstractVolume::createBoundarySurfaces()
{
  // transform Surfaces To BoundarySurfaces
  const std::vector<const Surface*> surfaces
      = Volume::volumeBounds().decomposeToSurfaces(m_transform);

  // counter to flip the inner/outer position for Cylinders
  int    sfCounter = 0;
  size_t sfNumber  = surfaces.size();

  for (auto& sf : surfaces) {
    // flip inner/outer for cylinders
    AbstractVolume* inner
        = (sf->type() == Surface::Cylinder && sfCounter == 3 && sfNumber > 3)
        ? nullptr
        : this;
    AbstractVolume* outer = (inner) ? nullptr : this;
    // create the boundary surface
    BoundarySurfacePtr bSurface
        = std::make_shared<const BoundarySurfaceT<AbstractVolume>>(
            std::unique_ptr<const Surface>(sf), inner, outer);
    m_boundarySurfaces.push_back(bSurface);
  }
}
