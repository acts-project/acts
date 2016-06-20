// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StrawSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/StrawSurface.hpp"
#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include <iomanip>
#include <iostream>

Acts::StrawSurface::StrawSurface(
    std::shared_ptr<Transform3D>    htrans,
    std::shared_ptr<const LineBounds> lbounds)
  : LineSurface(htrans,lbounds)
{
}

Acts::StrawSurface::StrawSurface(
    std::shared_ptr<const LineBounds> lbounds,
    const DetectorElementBase& detelement,
    const Identifier&                id)
  : LineSurface(lbounds,detelement, id)
{
}

Acts::StrawSurface::~StrawSurface()
{
}

Acts::StrawSurface&
Acts::StrawSurface::operator=(const StrawSurface& slsf)
{
  if (this != &slsf) {
    LineSurface::operator=(slsf);
    m_bounds  = slsf.m_bounds;
  }
  return *this;
}
