///////////////////////////////////////////////////////////////////
// PlaneSurface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// eometry module
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/RectangleBounds.h"
#include "Surfaces/TriangleBounds.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Surfaces/DiamondBounds.h"
#include "Surfaces/EllipseBounds.h"
#include "Surfaces/NoBounds.h"
// Core module
#include "Core/Identifier.h"
// STD/STL
#include <iostream>
#include <iomanip>

Acts::NoBounds Acts::PlaneSurface::s_boundless;

// default constructor
Acts::PlaneSurface::PlaneSurface() :
  Acts::Surface(),
  m_bounds(nullptr)
{}

// copy constructor
Acts::PlaneSurface::PlaneSurface(const PlaneSurface& psf) :
  Acts::Surface(psf),
  m_bounds(psf.m_bounds)
{}

// copy constructor with shift
Acts::PlaneSurface::PlaneSurface(const PlaneSurface& psf, const Acts::Transform3D& transf) :
  Acts::Surface(psf, transf),
  m_bounds(psf.m_bounds)
{}

// constructor from normal vector
Acts::PlaneSurface::PlaneSurface(const Acts::Vector3D& position, const Vector3D& normal) :
  Acts::Surface(),
  m_bounds(nullptr)
{
    Acts::Translation3D curvilinearTranslation(position.x(), position.y(),position.z());
    // the right-handed coordinate system is defined as
    // T = normal
    // U = Z x T if T not parallel to Z otherwise U = X x T
    // V = T x U
    // create the rotation
    Vector3D T = normal.normalized();
    Vector3D U = fabs(T.dot(Vector3D::UnitZ())) < 0.99 ? Vector3D::UnitZ().cross(T) : Vector3D::UnitX().cross(T);
    Vector3D V = T.cross(U);
    Acts::RotationMatrix3D curvilinearRotation;
    curvilinearRotation.col(0) = T;
    curvilinearRotation.col(1) = U;
    curvilinearRotation.col(2) = V;

    // curvilinear surfaces are boundless
    Acts::Surface::m_transform   = std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D);
    (*Acts::Surface::m_transform) = curvilinearRotation;
    Acts::Surface::m_transform->pretranslate(position);
}

// construct form DetectorElementBase & potentially identifier
Acts::PlaneSurface::PlaneSurface(const Acts::DetectorElementBase& detelement, const Identifier& identifier) :
  Acts::Surface(detelement, identifier),
  m_bounds(nullptr)
{}

// construct planar surface without bounds
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts:: Transform3D> htrans) :
  Acts::Surface(htrans),
  m_bounds(nullptr)
{}

// construct planar surface without bounds
Acts::PlaneSurface::PlaneSurface(std::unique_ptr<Acts::Transform3D> htrans) :
  Acts::Surface(std::move(htrans)),
  m_bounds(nullptr)
{}

// construct rectangle module
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts:: Transform3D> htrans, double halephi, double haleta) :
  Acts::Surface(htrans),
  m_bounds(new Acts::RectangleBounds(halephi, haleta))
{}

// construct trapezoidal module with parameters
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts:: Transform3D> htrans, double minhalephi, double maxhalephi, double haleta) :
  Acts::Surface(htrans),
  m_bounds(new Acts::TrapezoidBounds(minhalephi, maxhalephi, haleta))
{}

// construct module with shared boundaries
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts:: Transform3D> htrans, const Acts::PlanarBounds* tbounds) :
  Acts::Surface(htrans),
  m_bounds(tbounds)
{}

// construct module with shared boundaries
Acts::PlaneSurface::PlaneSurface(std::shared_ptr<Acts:: Transform3D> htrans, std::shared_ptr<const Acts::PlanarBounds> tbounds) :
  Acts::Surface(htrans),
  m_bounds(tbounds)
{}

// destructor (will call destructor from base class which deletes objects)
Acts::PlaneSurface::~PlaneSurface()
{}

Acts::PlaneSurface& Acts::PlaneSurface::operator=(const Acts::PlaneSurface& psf){

  if (this!=&psf){
    Acts::Surface::operator=(psf);
    m_bounds =  psf.m_bounds;
  }
  return *this;
}

bool Acts::PlaneSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::PlaneSurface* psf = dynamic_cast<const Acts::PlaneSurface*>(&sf);
  if (!psf) return false;
  if (psf==this) return true;
  bool transfEqual(transform().isApprox(psf->transform(), 10e-8));
  bool centerEqual = center() == psf->center();
  bool boundsEqual = bounds() == psf->bounds();
  return transfEqual&&centerEqual&&boundsEqual;
}

void Acts::PlaneSurface::localToGlobal(const Acts::Vector2D& locpos, const Acts::Vector3D&, Acts::Vector3D& glopos) const
{
	Acts::Vector3D loc3Dframe(locpos[Acts::eLOC_X], locpos[Acts::eLOC_Y], 0.);
    glopos = transform()*loc3Dframe;
}

bool Acts::PlaneSurface::globalToLocal(const Acts::Vector3D& glopos, const Acts::Vector3D&, Acts::Vector2D& locpos) const
{
	Acts::Vector3D loc3Dframe = (transform().inverse())*glopos;
    locpos = Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y());
    return (( loc3Dframe.z()*loc3Dframe.z() > s_onSurfaceTolerance*s_onSurfaceTolerance ) ? false : true );
}

bool Acts::PlaneSurface::isOnSurface(const Acts::Vector3D& glopo, const BoundaryCheck& bchk) const
{
	Acts::Vector3D loc3Dframe = (transform().inverse())*glopo;
    if ( fabs(loc3Dframe.z()) > s_onSurfaceTolerance ) return false;
    return ( bchk ?  bounds().inside(Acts::Vector2D(loc3Dframe.x(),loc3Dframe.y()),bchk) : true);
}
