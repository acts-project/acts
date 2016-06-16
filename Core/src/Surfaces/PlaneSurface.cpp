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
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/BoundlessT.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include <iomanip>
#include <iostream>

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& psf)
  : Surface(psf)
  , m_bounds(psf.m_bounds)
  , m_normal(psf.m_normal)
{}

Acts::PlaneSurface::PlaneSurface(const PlaneSurface& psf,
                                 const Transform3D&  transf)
  : Surface(psf, transf)
  , m_bounds(psf.m_bounds)
  , m_normal(0.,0.,0.) 
{
    m_normal = std::move(transform().rotation().col(2));
}

Acts::PlaneSurface::PlaneSurface(const Vector3D&  center,
                                 const Vector3D&  normal)
  : Surface()
  , m_bounds(nullptr)
  , m_normal(normal)
{
  Acts::Translation3D curvilinearTranslation(
      center.x(), center.y(), center.z());
  /// the right-handed coordinate system is defined as
  /// T = normal
  /// U = Z x T if T not parallel to Z otherwise U = X x T
  /// V = T x U
  Vector3D T = normal.normalized();
  Vector3D U = fabs(T.dot(Vector3D::UnitZ())) < 0.99
      ? Vector3D::UnitZ().cross(T)
      : Vector3D::UnitX().cross(T);
  Vector3D               V = T.cross(U);
  RotationMatrix3D       curvilinearRotation;
  curvilinearRotation.col(0) = T;
  curvilinearRotation.col(1) = U;
  curvilinearRotation.col(2) = V;

  // curvilinear surfaces are boundless
  Surface::m_transform    = std::make_shared<Transform3D>();
  (*Surface::m_transform) = curvilinearRotation;
  Surface::m_transform->pretranslate(center);
}

Acts::PlaneSurface::PlaneSurface( std::shared_ptr<const PlanarBounds> pbounds,
                                 const Acts::DetectorElementBase& detelement,
                                 const Identifier&                identifier)
  : Surface(detelement, identifier)
  , m_bounds(pbounds)
  , m_normal(0.,0.,0.)
{
    /// surfaces representing a detector element must have bounds
    assert(pbounds);
    m_normal = std::move(transform().rotation().col(2));    
}

Acts::PlaneSurface::PlaneSurface(
    std::shared_ptr<Transform3D>        htrans,
    std::shared_ptr<const PlanarBounds> pbounds)
  : Surface(std::move(htrans))
  , m_bounds(std::move(pbounds))
  , m_normal(0.,0.,0.) 
{
    m_normal = std::move(transform().rotation().col(2));    
}

Acts::PlaneSurface::~PlaneSurface()
{}

Acts::PlaneSurface&
Acts::PlaneSurface::operator=(const PlaneSurface& psf)
{
  if (this != &psf) {
    Surface::operator=(psf);
    m_bounds               = psf.m_bounds;
  }
  return *this;
}

void
Acts::PlaneSurface::localToGlobal(const Vector2D& lpos,
                                  const Vector3D&,
                                  Vector3D& gpos) const
{
  Acts::Vector3D loc3Dframe(lpos[Acts::eLOC_X], lpos[Acts::eLOC_Y], 0.);
  gpos = transform() * loc3Dframe;
}

bool
Acts::PlaneSurface::globalToLocal(const Vector3D& gpos,
                                  const Vector3D&,
                                  Acts::Vector2D& lpos) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * gpos;
  lpos                    = Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y());
  return ((loc3Dframe.z() * loc3Dframe.z()
           > s_onSurfaceTolerance * s_onSurfaceTolerance)
              ? false
              : true);
}

bool
Acts::PlaneSurface::isOnSurface(const Vector3D& glopo,
                                const BoundaryCheck&  bchk) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * glopo;
  if (fabs(loc3Dframe.z()) > s_onSurfaceTolerance) return false;
  return (bchk
              ? bounds().inside(Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y()),
                                bchk)
              : true);
}
