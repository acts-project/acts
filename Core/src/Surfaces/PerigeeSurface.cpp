// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// PerigeeSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/PerigeeSurface.hpp"

#include <iomanip>
#include <iostream>

Acts::PerigeeSurface::PerigeeSurface(const Vector3D& gp) : LineSurface(nullptr)
{
  Surface::m_transform = std::make_shared<const Transform3D>(
      Translation3D(gp.x(), gp.y(), gp.z()));
}

Acts::PerigeeSurface::PerigeeSurface(
    std::shared_ptr<const Transform3D> tTransform)
  : LineSurface(tTransform)
{
}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& other)
  : LineSurface(other)
{
}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& other,
                                     const Transform3D&    shift)
  : LineSurface(other, shift)
{
}

Acts::PerigeeSurface::~PerigeeSurface()
{
}

Acts::PerigeeSurface&
Acts::PerigeeSurface::operator=(const PerigeeSurface& other)
{
  if (this != &other) {
    LineSurface::operator=(other);
  }
  return *this;
}

Acts::PerigeeSurface*
Acts::PerigeeSurface::clone(const Transform3D* shift) const
{
  if (shift) return new PerigeeSurface(*this, *shift);
  return new PerigeeSurface(*this);
}

Acts::Surface::SurfaceType
Acts::PerigeeSurface::type() const
{
  return Surface::Perigee;
}

std::string
Acts::PerigeeSurface::name() const
{
  return "Acts::PerigeeSurface";
}

std::ostream&
Acts::PerigeeSurface::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::PerigeeSurface:" << std::endl;
  sl << "     Center position  (x, y, z) = (" << center().x() << ", "
     << center().y() << ", " << center().z() << ")";
  sl << std::setprecision(-1);
  return sl;
}
