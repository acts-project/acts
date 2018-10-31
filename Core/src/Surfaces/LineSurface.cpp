// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LineSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/LineSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               double                             radius,
                               double                             halez)
  : GeometryObject()
  , Surface(std::move(htrans))
  , m_bounds(std::make_shared<const LineBounds>(radius, halez))
{
}

Acts::LineSurface::LineSurface(std::shared_ptr<const Transform3D> htrans,
                               std::shared_ptr<const LineBounds>  lbounds)
  : GeometryObject(), Surface(std::move(htrans)), m_bounds(std::move(lbounds))
{
}

Acts::LineSurface::LineSurface(const std::shared_ptr<const LineBounds>& lbounds,
                               const DetectorElementBase& detelement)
  : GeometryObject(), Surface(detelement), m_bounds(lbounds)
{
  throw_assert(lbounds, "LineBounds must not be nullptr");
}

Acts::LineSurface::LineSurface(const LineSurface& other)
  : GeometryObject(), Surface(other), m_bounds(other.m_bounds)
{
}

Acts::LineSurface::LineSurface(const LineSurface& other,
                               const Transform3D& transf)
  : GeometryObject(), Surface(other, transf), m_bounds(other.m_bounds)
{
}

Acts::LineSurface&
Acts::LineSurface::operator=(const LineSurface& other)
{
  if (this != &other) {
    Surface::operator=(other);
    m_bounds         = other.m_bounds;
  }
  return *this;
}

Acts::LineSurface::LineSurface(const variant_data& vardata) : GeometryObject()
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "LineSurface", "Variant data type must be LineSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map bounds     = payload.get<variant_map>("bounds");
  std::string boundsType = bounds.get<std::string>("type");

  throw_assert(boundsType == "LineBounds",
               "Can only construct LineSurface from LineBounds");

  m_bounds = std::make_shared<const LineBounds>(bounds);

  if (payload.count("transform") != 0u) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
}

Acts::variant_data
Acts::LineSurface::toVariantData() const
{
  using namespace std::string_literals;

  variant_map payload;

  variant_data bounds = m_bounds->toVariantData();
  payload["bounds"]   = bounds;

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  variant_map data;
  data["type"]    = "LineSurface"s;
  data["payload"] = payload;
  return data;
}
