// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedDiscSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/SubtractedDiscSurface.hpp"

Acts::SubtractedDiscSurface::SubtractedDiscSurface(
    const SubtractedDiscSurface& psf)
  : Acts::DiscSurface(psf), m_subtrVol(psf.m_subtrVol), m_shared(psf.m_shared)
{
}

Acts::SubtractedDiscSurface::SubtractedDiscSurface(
    const SubtractedDiscSurface& psf,
    const Acts::Transform3D&     shift)
  : Acts::DiscSurface(psf, shift)
  , m_subtrVol(psf.m_subtrVol)
  , m_shared(psf.m_shared)
{
}

Acts::SubtractedDiscSurface::SubtractedDiscSurface(const Acts::DiscSurface& ps,
                                                   AreaExcluder*            vol,
                                                   bool shared)
  : Acts::DiscSurface(ps)
  , m_subtrVol(std::shared_ptr<AreaExcluder>(vol))
  , m_shared(shared)
{
}

Acts::SubtractedDiscSurface::~SubtractedDiscSurface()
{
}

Acts::SubtractedDiscSurface&
Acts::SubtractedDiscSurface::operator=(const Acts::SubtractedDiscSurface& psf)
{
  if (this != &psf) {
    Acts::DiscSurface::operator=(psf);
    m_subtrVol                 = psf.m_subtrVol;
    m_shared                   = psf.m_shared;
  }
  return *this;
}

bool
Acts::SubtractedDiscSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::SubtractedDiscSurface* sdsf
      = dynamic_cast<const Acts::SubtractedDiscSurface*>(&sf);
  if (!sdsf) return false;
  bool surfaceEqual = Acts::DiscSurface::operator==(sf);
  bool sharedEqual  = (surfaceEqual) ? (shared() == sdsf->shared()) : false;
  return sharedEqual;
}
