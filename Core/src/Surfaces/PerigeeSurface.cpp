// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// PerigeeSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

Acts::PerigeeSurface::PerigeeSurface(const Vector3D& gp)
  : LineSurface(nullptr, nullptr)
{
  Surface::m_transform = std::make_shared<const Transform3D>(
      Translation3D(gp.x(), gp.y(), gp.z()));
}

Acts::PerigeeSurface::PerigeeSurface(
    std::shared_ptr<const Transform3D> tTransform)
  : GeometryObject(), LineSurface(std::move(tTransform))
{
}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& other)
  : GeometryObject(), LineSurface(other)
{
}

}

Acts::PerigeeSurface&
Acts::PerigeeSurface::operator=(const PerigeeSurface& other)
{
  if (this != &other) {
    LineSurface::operator=(other);
  }
  return *this;
}

std::shared_ptr<Acts::PerigeeSurface>
Acts::PerigeeSurface::clone() const
{
  return std::shared_ptr<PerigeeSurface>(this->clone_impl());
}

Acts::PerigeeSurface*
Acts::PerigeeSurface::clone_impl() const
{
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
Acts::PerigeeSurface::toStream(Context ctx, std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::PerigeeSurface:" << std::endl;
  const Vector3D& sfCenter = center(ctx);
  sl << "     Center position  (x, y, z) = (" << sfCenter.x() << ", "
     << sfCenter.y() << ", " << sfCenter.z() << ")";
  sl << std::setprecision(-1);
  return sl;
}
