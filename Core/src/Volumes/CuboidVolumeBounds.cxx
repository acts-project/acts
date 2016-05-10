///////////////////////////////////////////////////////////////////
// CuboidVolumeBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/CuboidVolumeBounds.h"
#include "ACTS/Surfaces/Surface.h"
#include "ACTS/Surfaces/PlaneSurface.h"
#include "ACTS/Surfaces/RectangleBounds.h"
// STD/STL
#include <iostream>
#include <math.h>

Acts::CuboidVolumeBounds::CuboidVolumeBounds() :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{}

Acts::CuboidVolumeBounds::CuboidVolumeBounds(double halex, double haley, double halez) :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{
  m_boundValues[bv_halfX] = halex;
  m_boundValues[bv_halfY] = haley;
  m_boundValues[bv_halfZ] = halez;
}

Acts::CuboidVolumeBounds::CuboidVolumeBounds(const Acts::CuboidVolumeBounds& bobo) :
 VolumeBounds(),
 m_boundValues(bobo.m_boundValues)
{}

Acts::CuboidVolumeBounds::~CuboidVolumeBounds()
{}

Acts::CuboidVolumeBounds& Acts::CuboidVolumeBounds::operator=(const Acts::CuboidVolumeBounds& bobo)
{
  if (this!=&bobo)
      m_boundValues = bobo.m_boundValues;
  return *this;
}

const std::vector<const Acts::Surface*>* Acts::CuboidVolumeBounds::decomposeToSurfaces(std::shared_ptr<Acts::Transform3D> transformPtr) const
{

    // the transform
    Acts::Transform3D transform = ( transformPtr == nullptr) ? Acts::Transform3D::Identity() : (*transformPtr.get());
    Acts::Transform3D* tTransform = 0;

    std::vector<const Acts::Surface*>* retsf = new std::vector<const Acts::Surface*>;
    // memory optimisation
    retsf->reserve(6);
    // face surfaces xy -------------------------------------
    //   (1) - at negative local z
    std::shared_ptr<const Acts::PlanarBounds> xyBounds(faceXYRectangleBounds());
    tTransform = new Acts::Transform3D( transform * Acts::AngleAxis3D(M_PI, Acts::Vector3D(0.,1.,0.))*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), xyBounds));
    //   (2) - at positive local z
    tTransform =  new Acts::Transform3D( transform*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), xyBounds));
    // face surfaces yz -------------------------------------
    // transmute cyclical
    //   (3) - at negative local x
    std::shared_ptr<const Acts::PlanarBounds> yzBounds(faceYZRectangleBounds());
    tTransform =  new Acts::Transform3D( transform  *Acts::AngleAxis3D(M_PI, Acts::Vector3D(0.,0.,1.))
								                   *Acts::Translation3D(Acts::Vector3D(halflengthX(), 0.,0))
								                   *Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0.,1.,0))*Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0.,0.,1.)));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), yzBounds));
    //   (4) - at positive local x
    tTransform = new Acts::Transform3D( transform *Acts::Translation3D(Acts::Vector3D(halflengthX(),0.,0.))
								                *Acts::AngleAxis3D(0.5 * M_PI,Acts::Vector3D(0.,1.,0.))*Acts::AngleAxis3D(0.5 * M_PI, Acts::Vector3D(0.,0.,1.)));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), yzBounds));
    // face surfaces zx -------------------------------------
    std::shared_ptr<const Acts::PlanarBounds> zxBounds(faceZXRectangleBounds());
    //   (5) - at negative local y
    tTransform = new Acts::Transform3D( transform *Acts::AngleAxis3D(M_PI, Acts::Vector3D(1.,0.,0.))
								                 *Acts::Translation3D(Acts::Vector3D(0.,halflengthY(),0.))
								                 *Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0.,1.,0.))*Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(1.,0.,0.)));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), zxBounds));
    //   (6) - at positive local y
    tTransform = new Acts::Transform3D( transform*Acts::Translation3D(Acts::Vector3D(0.,halflengthY(),0.))
								                *Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(0.,1.,0.))*Acts::AngleAxis3D(-0.5 * M_PI, Acts::Vector3D(1.,0.,0.)));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), zxBounds));
    // return the surfaces
    return retsf;
}

Acts::RectangleBounds* Acts::CuboidVolumeBounds::faceXYRectangleBounds() const
{ return new Acts::RectangleBounds(m_boundValues[bv_halfX], m_boundValues[bv_halfY]); }

Acts::RectangleBounds* Acts::CuboidVolumeBounds::faceYZRectangleBounds() const
{ return new Acts::RectangleBounds(m_boundValues[bv_halfY], m_boundValues[bv_halfZ]); }

Acts::RectangleBounds* Acts::CuboidVolumeBounds::faceZXRectangleBounds() const
{ return new Acts::RectangleBounds(m_boundValues[bv_halfZ], m_boundValues[bv_halfX]); }

// ostream operator overload
std::ostream& Acts::CuboidVolumeBounds::dump( std::ostream& sl ) const
{
    return dumpT(sl);
}
