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

// default constructor
Acts::PlaneSurface::PlaneSurface() 
    : Surface()
    , m_bounds(nullptr)
    , m_normal(0.,0.,0.)    

{}

// copy constructor
Acts::PlaneSurface::PlaneSurface(const PlaneSurface& psf)
  : Surface(psf)
  , m_bounds(psf.m_bounds)
  , m_normal(psf.m_normal)
{}

// copy constructor with shift
Acts::PlaneSurface::PlaneSurface(const PlaneSurface&      psf,
                                 const Acts::Transform3D& transf)
  : Surface(psf, transf)
  , m_bounds(psf.m_bounds)
  , m_normal(0.,0.,0.) 
{
    m_normal = transform().rotation().col(2);
}

// constructor from normal vector
Acts::PlaneSurface::PlaneSurface(const Acts::Vector3D& position,
                                 const Vector3D&       normal)
  : Surface()
  , m_bounds(nullptr)
  , m_normal(normal)
{
  Acts::Translation3D curvilinearTranslation(
      position.x(), position.y(), position.z());
  // the right-handed coordinate system is defined as
  // T = normal
  // U = Z x T if T not parallel to Z otherwise U = X x T
  // V = T x U
  // create the rotation
  Vector3D T = normal.normalized();
  Vector3D U = fabs(T.dot(Vector3D::UnitZ())) < 0.99
      ? Vector3D::UnitZ().cross(T)
      : Vector3D::UnitX().cross(T);
  Vector3D               V = T.cross(U);
  Acts::RotationMatrix3D curvilinearRotation;
  curvilinearRotation.col(0) = T;
  curvilinearRotation.col(1) = U;
  curvilinearRotation.col(2) = V;

  // curvilinear surfaces are boundless
  Surface::m_transform    = std::make_shared<Acts::Transform3D>();
  (*Surface::m_transform) = curvilinearRotation;
  Surface::m_transform->pretranslate(position);
}

// construct form DetectorElementBase & potentially identifier
Acts::PlaneSurface::PlaneSurface(const Acts::DetectorElementBase& detelement,
                                 const Identifier&                identifier)
  : Surface(detelement, identifier)
  , m_bounds(nullptr)
  , m_normal(0.,0.,0.)
{
    m_normal = transform().rotation().col(2);    
}

// construct planar surface without bounds
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts::Transform3D> htrans)
  : Surface(htrans)
  , m_bounds(std::make_shared<BoundlessT<RectangleBounds>()>)
  , m_normal(0.,0.,0.) 
{
    m_normal = transform().rotation().col(2);    
}

// construct planar surface without bounds
Acts::PlaneSurface::PlaneSurface(std::unique_ptr<Acts::Transform3D> htrans)
  : Surface(std::move(htrans))
  , m_bounds(nullptr)
  , m_normal(0.,0.,0.) 
{
    m_normal = transform().rotation().col(2);    
}

// construct module with shared boundaries
Acts::PlaneSurface::PlaneSurface(
    std::shared_ptr<Acts::Transform3D>        htrans,
    std::shared_ptr<const Acts::PlanarBounds> tbounds)
  : Surface(std::move(htrans))
  , m_bounds(nullptr)
  , m_normal(0.,0.,0.) 
{
    m_normal = transform().rotation().col(2);    
}

// destructor (will call destructor from base class which deletes objects)
Acts::PlaneSurface::~PlaneSurface()
{}

Acts::PlaneSurface&
Acts::PlaneSurface::operator=(const Acts::PlaneSurface& psf)
{
  if (this != &psf) {
    Surface::operator=(psf);
    m_bounds               = psf.m_bounds;
  }
  return *this;
}

bool
Acts::PlaneSurface::operator==(const Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const PlaneSurface* psf = dynamic_cast<const PlaneSurface*>(&sf);
  if (!psf) return false;
  if (psf == this) return true;
  // @TODO make approx parameter more useful
  bool transfEqual(transform().isApprox(psf->transform(), 10e-8));
  bool centerEqual = center() == psf->center();
  bool boundsEqual = bounds() == psf->bounds();
  return transfEqual && centerEqual && boundsEqual;
}

void
Acts::PlaneSurface::localToGlobal(const Acts::Vector2D& locpos,
                                  const Acts::Vector3D&,
                                  Acts::Vector3D& glopos) const
{
  Acts::Vector3D loc3Dframe(locpos[Acts::eLOC_X], locpos[Acts::eLOC_Y], 0.);
  glopos = transform() * loc3Dframe;
}

bool
Acts::PlaneSurface::globalToLocal(const Acts::Vector3D& glopos,
                                  const Acts::Vector3D&,
                                  Acts::Vector2D& locpos) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * glopos;
  locpos                    = Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y());
  return ((loc3Dframe.z() * loc3Dframe.z()
           > s_onSurfaceTolerance * s_onSurfaceTolerance)
              ? false
              : true);
}

bool
Acts::PlaneSurface::isOnSurface(const Acts::Vector3D& glopo,
                                const BoundaryCheck&  bchk) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse()) * glopo;
  if (fabs(loc3Dframe.z()) > s_onSurfaceTolerance) return false;
  return (bchk
              ? bounds().inside(Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y()),
                                bchk)
              : true);
}
