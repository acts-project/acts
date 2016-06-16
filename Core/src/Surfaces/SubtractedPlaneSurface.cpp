// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedPlaneSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/SubtractedPlaneSurface.hpp"

Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(
    const SubtractedPlaneSurface& psf)
  : Acts::PlaneSurface(psf), m_subtrVol(psf.m_subtrVol), m_shared(psf.m_shared)
{
}

Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(
    const SubtractedPlaneSurface& psf,
    const Acts::Transform3D&      transf)
  : Acts::PlaneSurface(psf, transf)
  , m_subtrVol(psf.m_subtrVol)
  , m_shared(psf.m_shared)
{
}

Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(
    const Acts::PlaneSurface& ps,
    AreaExcluder*             vol,
    bool                      shared)
  : Acts::PlaneSurface(ps), m_subtrVol(vol), m_shared(shared)
{
}

Acts::SubtractedPlaneSurface::~SubtractedPlaneSurface()
{
}

Acts::SubtractedPlaneSurface&
Acts::SubtractedPlaneSurface::operator=(const SubtractedPlaneSurface& psf)
{
  if (this != &psf) {
    Acts::PlaneSurface::operator=(psf);
    m_subtrVol                  = psf.m_subtrVol;
    m_shared                    = psf.m_shared;
  }
  return *this;
}

bool
Acts::SubtractedPlaneSurface::operator==(const Acts::Surface& sf) const
{
  // fast exit
  if (&sf == this) return true;
  // first check the type not to compare apples with oranges
  const Acts::SubtractedPlaneSurface* spsf
      = dynamic_cast<const Acts::SubtractedPlaneSurface*>(&sf);
  if (!spsf) return false;
  bool surfaceEqual = Acts::PlaneSurface::operator==(sf);
  bool sharedEqual  = (surfaceEqual) ? (shared() == spsf->shared()) : false;
  return sharedEqual;
}
