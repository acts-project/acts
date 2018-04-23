// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StrawSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/StrawSurface.hpp"

#include <iomanip>
#include <iostream>

#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/VariantData.hpp"

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 double                             radius,
                                 double                             halez)
  : GeometryObject(), LineSurface(htrans, radius, halez)
{
}

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 std::shared_ptr<const LineBounds>  lbounds)
  : GeometryObject(), LineSurface(htrans, lbounds)
{
}

Acts::StrawSurface::StrawSurface(std::shared_ptr<const LineBounds> lbounds,
                                 const DetectorElementBase&        detelement,
                                 const Identifier&                 id)
  : GeometryObject(), LineSurface(lbounds, detelement, id)
{
}

Acts::StrawSurface::StrawSurface(const Acts::StrawSurface& other)
  : GeometryObject(), LineSurface(other)
{
}

Acts::StrawSurface::StrawSurface(const StrawSurface& other,
                                 const Transform3D&  htrans)
  : LineSurface(other, htrans)
{
}

Acts::StrawSurface::StrawSurface(const variant_data& data_)
  : GeometryObject(), LineSurface(nullptr, nullptr)
{
  throw_assert(data_.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(data_);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "StrawSurface",
               "Variant data type must be StrawSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map bounds     = payload.get<variant_map>("bounds");
  std::string boundsType = bounds.get<std::string>("type");

  throw_assert(boundsType == "LineBounds",
               "Can only construct StrawSurface from LineBounds");

  m_bounds = std::make_shared<const LineBounds>(bounds);

  if (payload.count("transform")) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
}

Acts::StrawSurface::~StrawSurface()
{
}

Acts::StrawSurface&
Acts::StrawSurface::operator=(const StrawSurface& other)
{
  if (this != &other) {
    LineSurface::operator=(other);
    m_bounds             = other.m_bounds;
  }
  return *this;
}

Acts::StrawSurface*
Acts::StrawSurface::clone(const Transform3D* shift) const
{
  if (shift) new StrawSurface(*this, *shift);
  return new StrawSurface(*this);
}

Acts::variant_data
Acts::StrawSurface::toVariantData() const
{
  using namespace std::string_literals;

  variant_map  payload;
  variant_data bounds = m_bounds->toVariantData();
  payload["bounds"]   = bounds;
  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  variant_map data;
  data["type"]    = "StrawSurface"s;
  data["payload"] = payload;
  return data;
}
