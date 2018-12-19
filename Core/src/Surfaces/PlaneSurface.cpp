// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/PlaneSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include "Acts/Utilities/InstanceFactory.hpp"

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other)
  : GeometryObject(), Surface(other), m_bounds(other.m_bounds)
{
}

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other,
                                 const Transform3D&  transf)
  : GeometryObject(), Surface(other, transf), m_bounds(other.m_bounds)
{
}

Acts::PlaneSurface::PlaneSurface(const Vector3D& center, const Vector3D& normal)
  : Surface(), m_bounds(nullptr)
{
  /// the right-handed coordinate system is defined as
  /// T = normal
  /// U = Z x T if T not parallel to Z otherwise U = X x T
  /// V = T x U
  Vector3D T = normal.normalized();
  Vector3D U = std::abs(T.dot(Vector3D::UnitZ())) < s_curvilinearProjTolerance
      ? Vector3D::UnitZ().cross(T).normalized()
      : Vector3D::UnitX().cross(T).normalized();
  Vector3D         V = T.cross(U);
  RotationMatrix3D curvilinearRotation;
  curvilinearRotation.col(0) = U;
  curvilinearRotation.col(1) = V;
  curvilinearRotation.col(2) = T;

  // curvilinear surfaces are boundless
  Transform3D transform{curvilinearRotation};
  transform.pretranslate(center);
  Surface::m_transform = std::make_shared<const Transform3D>(transform);
}

Acts::PlaneSurface::PlaneSurface(
    const std::shared_ptr<const PlanarBounds>& pbounds,
    const Acts::DetectorElementBase&           detelement)
  : Surface(detelement), m_bounds(pbounds)
{
  /// surfaces representing a detector element must have bounds
  throw_assert(pbounds, "PlaneBounds must not be nullptr");
}

Acts::PlaneSurface::PlaneSurface(std::shared_ptr<const Transform3D>  htrans,
                                 std::shared_ptr<const PlanarBounds> pbounds)
  : Surface(std::move(htrans)), m_bounds(std::move(pbounds))
{
}

Acts::PlaneSurface::PlaneSurface(const variant_data& vardata)
{
  // we need to figure out which way the PS was constructed before
  throw_assert(vardata.which() == 4, "Variant data must be map");
  variant_map data = boost::get<variant_map>(vardata);
  throw_assert(data.count("type"), "Variant data must have type.");
  // std::string type = boost::get<std::string>(data["type"]);
  std::string type = data.get<std::string>("type");
  throw_assert(type == "PlaneSurface",
               "Variant data type must be PlaneSurface");

  variant_map payload    = data.get<variant_map>("payload");
  variant_map bounds     = payload.get<variant_map>("bounds");
  std::string boundsType = bounds.get<std::string>("type");

  InstanceFactory                     factory;
  std::shared_ptr<const PlanarBounds> pbounds
      = factory.planarBounds(boundsType, bounds);

  m_bounds = pbounds;

  if (payload.count("transform") != 0u) {
    // we have a transform
    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));
    m_transform = trf;
  }
}

Acts::PlaneSurface&
Acts::PlaneSurface::operator=(const PlaneSurface& other)
{
  if (this != &other) {
    Surface::operator=(other);
    m_bounds         = other.m_bounds;
  }
  return *this;
}

Acts::Surface::SurfaceType
Acts::PlaneSurface::type() const
{
  return Surface::Plane;
}

void
Acts::PlaneSurface::localToGlobal(const Vector2D& lpos,
                                  const Vector3D& /*gmom*/,
                                  Vector3D& gpos) const
{
  Vector3D loc3Dframe(lpos[Acts::eLOC_X], lpos[Acts::eLOC_Y], 0.);
  /// the chance that there is no transform is almost 0, let's apply it
  gpos = transform() * loc3Dframe;
}

bool
Acts::PlaneSurface::globalToLocal(const Vector3D& gpos,
                                  const Vector3D& /*gmom*/,
                                  Acts::Vector2D& lpos) const
{
  /// the chance that there is no transform is almost 0, let's apply it
  Vector3D loc3Dframe = (transform().inverse()) * gpos;
  lpos                = Vector2D(loc3Dframe.x(), loc3Dframe.y());
  return ((loc3Dframe.z() * loc3Dframe.z()
           > s_onSurfaceTolerance * s_onSurfaceTolerance)
              ? false
              : true);
}

std::string
Acts::PlaneSurface::name() const
{
  return "Acts::PlaneSurface";
}

std::shared_ptr<Acts::PlaneSurface>
Acts::PlaneSurface::clone(const Transform3D* shift) const
{
  return std::shared_ptr<PlaneSurface>(this->clone_impl(shift));
}

Acts::PlaneSurface*
Acts::PlaneSurface::clone_impl(const Transform3D* shift) const
{
  if (shift != nullptr) {
    return new PlaneSurface(*this, *shift);
  }
  return new PlaneSurface(*this);
}

const Acts::SurfaceBounds&
Acts::PlaneSurface::bounds() const
{
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

Acts::variant_data
Acts::PlaneSurface::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;

  variant_data bounds = m_bounds->toVariantData();
  payload["bounds"]   = bounds;

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  } else if (m_associatedDetElement != nullptr) {
    payload["transform"] = to_variant(m_associatedDetElement->transform());
  }

  variant_map data;
  data["type"]    = "PlaneSurface"s;
  data["payload"] = payload;

  return data;
}
