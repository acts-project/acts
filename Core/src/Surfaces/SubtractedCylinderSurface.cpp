// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SubtractedCylinderSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/SubtractedCylinderSurface.hpp"

// default constructor
Acts::SubtractedCylinderSurface::SubtractedCylinderSurface() :
  Acts::CylinderSurface(),
  m_subtrVol(),
  m_shared(true)
{}

// copy constructor
Acts::SubtractedCylinderSurface::SubtractedCylinderSurface(const SubtractedCylinderSurface& psf) :
  Acts::CylinderSurface(psf),
  m_subtrVol(psf.m_subtrVol),
  m_shared(psf.m_shared)
{}

// copy constructor with shift
Acts::SubtractedCylinderSurface::SubtractedCylinderSurface(const SubtractedCylinderSurface& psf, const Acts::Transform3D& transf) :
  Acts::CylinderSurface(psf, transf),
  m_subtrVol(psf.m_subtrVol),
  m_shared(psf.m_shared)
{}

// constructor
Acts::SubtractedCylinderSurface::SubtractedCylinderSurface(const Acts::CylinderSurface& ps, AreaExcluder* vol, bool shared) :
  Acts::CylinderSurface(ps),
  m_subtrVol(std::shared_ptr<Acts::AreaExcluder>(vol)),
  m_shared(shared)
{}

// destructor (will call destructor from base class which deletes objects)
Acts::SubtractedCylinderSurface::~SubtractedCylinderSurface()
{}

Acts::SubtractedCylinderSurface& Acts::SubtractedCylinderSurface::operator=(const Acts::SubtractedCylinderSurface& psf){

  if (this!=&psf){
    Acts::CylinderSurface::operator=(psf);
    m_subtrVol = psf.m_subtrVol;
    m_shared   = psf.m_shared;
  }
  return *this;

}

bool Acts::SubtractedCylinderSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::SubtractedCylinderSurface* scsf = dynamic_cast<const Acts::SubtractedCylinderSurface*>(&sf);
  if (!scsf) return false;
  bool surfaceEqual = Acts::CylinderSurface::operator==(sf);
  bool sharedEqual  = (surfaceEqual) ? (shared() == scsf->shared()) : false;
  return sharedEqual;
}

