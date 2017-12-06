// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/PlaneSurface.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Identifier.hpp"

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other)
  : Surface(other), m_bounds(other.m_bounds)
{
}

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& other,
                                 const Transform3D&  transf)
  : Surface(other, transf), m_bounds(other.m_bounds)
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

Acts::PlaneSurface::PlaneSurface(std::shared_ptr<const PlanarBounds> pbounds,
                                 const Acts::DetectorElementBase&    detelement,
                                 const Identifier&                   identifier)
  : Surface(detelement, identifier), m_bounds(pbounds)
{
  /// surfaces representing a detector element must have bounds
  assert(pbounds);
}

Acts::PlaneSurface::PlaneSurface(std::shared_ptr<const Transform3D>  htrans,
                                 std::shared_ptr<const PlanarBounds> pbounds)
  : Surface(std::move(htrans)), m_bounds(std::move(pbounds))
{
}

Acts::PlaneSurface::~PlaneSurface()
{
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
                                  const Vector3D&,
                                  Vector3D& gpos) const
{
  Vector3D loc3Dframe(lpos[Acts::eLOC_X], lpos[Acts::eLOC_Y], 0.);
  /// the chance that there is no transform is almost 0, let's apply it
  gpos = transform() * loc3Dframe;
}

bool
Acts::PlaneSurface::globalToLocal(const Vector3D& gpos,
                                  const Vector3D&,
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

bool
Acts::PlaneSurface::isOnSurface(const Vector3D&      glopo,
                                const BoundaryCheck& bcheck) const
{
  /// the chance that there is no transform is almost 0, let's apply it
  Vector3D loc3Dframe = (transform().inverse()) * glopo;
  if (std::abs(loc3Dframe.z()) > s_onSurfaceTolerance) return false;
  return (
      bcheck ? bounds().inside(Vector2D(loc3Dframe.x(), loc3Dframe.y()), bcheck)
             : true);
}

Acts::PlaneSurface*
Acts::PlaneSurface::clone(const Transform3D* shift) const
{
  if (shift) new PlaneSurface(*this, *shift);
  return new PlaneSurface(*this);
}

const Acts::SurfaceBounds&
Acts::PlaneSurface::bounds() const
{
  if (m_bounds) return (*m_bounds.get());
  return s_noBounds;
}

const Acts::Vector3D
Acts::PlaneSurface::normal(const Acts::Vector2D&) const
{
  // fast access via tranform matrix (and not rotation())
  auto tMatrix = transform().matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

const Acts::Vector3D
    Acts::PlaneSurface::binningPosition(Acts::BinningValue) const
{
  return center();
}

double
Acts::PlaneSurface::pathCorrection(const Acts::Vector3D&,
                                   const Acts::Vector3D& mom) const
{
  /// we can ignore the global position here
  return 1. / std::abs(normal().dot(mom.unit()));
}

Acts::Intersection
Acts::PlaneSurface::intersectionEstimate(
    const Acts::Vector3D&      gpos,
    const Acts::Vector3D&      dir,
    bool                       forceDir,
    const Acts::BoundaryCheck& bcheck) const
{
  double denom = dir.dot(normal());
  if (denom) {
    double   u = (normal().dot((center() - gpos))) / (denom);
    Vector3D intersectPoint(gpos + u * dir);
    // evaluate the intersection in terms of direction
    bool isValid = forceDir ? (u > 0.) : true;
    // evaluate (if necessary in terms of boundaries)
    isValid
        = bcheck ? (isValid && isOnSurface(intersectPoint, bcheck)) : isValid;
    // return the result
    return Intersection(intersectPoint, u, isValid);
  }
  return Intersection(gpos, 0., false);
}
