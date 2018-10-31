// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StrawSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Utilities/VariantData.hpp"

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 double                             radius,
                                 double                             halez)
  : GeometryObject(), LineSurface(std::move(htrans), radius, halez)
{
}

Acts::StrawSurface::StrawSurface(std::shared_ptr<const Transform3D> htrans,
                                 std::shared_ptr<const LineBounds>  lbounds)
  : GeometryObject(), LineSurface(std::move(htrans), std::move(lbounds))
{
}

Acts::StrawSurface::StrawSurface(
    const std::shared_ptr<const LineBounds>& lbounds,
    const DetectorElementBase&               detelement)
  : GeometryObject(), LineSurface(lbounds, detelement)
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

Acts::StrawSurface::StrawSurface(const variant_data& vardata)
  : GeometryObject(), LineSurface(nullptr, nullptr)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
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

  if (payload.count("transform") != 0u) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
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

std::shared_ptr<Acts::StrawSurface>
Acts::StrawSurface::clone(const Transform3D* shift) const
{
  return std::shared_ptr<StrawSurface>(this->clone_impl(shift));
}

Acts::StrawSurface*
Acts::StrawSurface::clone_impl(const Transform3D* shift) const
{
  if (shift != nullptr) {
    return new StrawSurface(*this, *shift);
  }
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

Acts::PolyhedronRepresentation
Acts::StrawSurface::polyhedronRepresentation(size_t l0div,
                                             size_t /*l1div*/) const
{
  std::vector<Vector3D>            vertices;
  std::vector<std::vector<size_t>> faces;

  if (l0div == 1) {
    throw std::domain_error("Polyhedron repr of straw with 1 div is undefined");
  }

  double phistep = 2 * M_PI / l0div;
  double hlZ     = m_bounds->halflengthZ();
  double r       = m_bounds->r();

  Vector3D left(r, 0, -hlZ);
  Vector3D right(r, 0, hlZ);

  for (size_t i = 0; i < l0div; i++) {
    Transform3D rot(AngleAxis3D(i * phistep, Vector3D::UnitZ()));
    vertices.push_back(transform() * rot * left);
    vertices.push_back(transform() * rot * right);
  }

  for (size_t v = 0; v < vertices.size() - 2; v = v + 2) {

    faces.push_back({v, v + 1, v + 3, v + 2});
  }
  if (l0div > 2) {
    faces.push_back({vertices.size() - 2, vertices.size() - 1, 1, 0});
  }

  return PolyhedronRepresentation(vertices, faces);
}
