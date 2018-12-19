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
#include "Acts/Utilities/VariantData.hpp"

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

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& other,
                                     const Transform3D&    shift)
  : GeometryObject(), LineSurface(other, shift)
{
}

Acts::PerigeeSurface::PerigeeSurface(const variant_data& vardata)
  : GeometryObject(), LineSurface(nullptr, nullptr)
{

  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "PerigeeSurface",
               "Variant data type must be PerigeeSurface");

  variant_map payload = data.get<variant_map>("payload");

  if (payload.count("transform") != 0u) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
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
Acts::PerigeeSurface::clone(const Transform3D* shift) const
{
  return std::shared_ptr<PerigeeSurface>(this->clone_impl(shift));
}

Acts::PerigeeSurface*
Acts::PerigeeSurface::clone_impl(const Transform3D* shift) const
{
  if (shift != nullptr) {
    return new PerigeeSurface(*this, *shift);
  }
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

Acts::variant_data
Acts::PerigeeSurface::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  variant_map data;
  data["type"]    = "PerigeeSurface"s;
  data["payload"] = payload;
  return data;
}
