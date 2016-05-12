///////////////////////////////////////////////////////////////////
// CylinderVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
// STD/STL
#include <iostream>
#include <math.h>

double Acts::CylinderVolumeBounds::s_numericalStable = 10e-2;

Acts::CylinderVolumeBounds::CylinderVolumeBounds() :
 VolumeBounds(),
 m_boundValues(4,0.)
{}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double radius, double halez) :
 VolumeBounds(),
 m_boundValues(4,0.)
{
    m_boundValues[bv_innerRadius]   = 0.;
    m_boundValues[bv_outerRadius]   = fabs(radius);
    m_boundValues[bv_halfPhiSector] = M_PI;
    m_boundValues[bv_halfZ]         = fabs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner, double router, double halez) :
 VolumeBounds(),
 m_boundValues(4,0.)
{
    m_boundValues[bv_innerRadius]   = fabs(rinner);
    m_boundValues[bv_outerRadius]   = fabs(router);
    m_boundValues[bv_halfPhiSector] = M_PI;
    m_boundValues[bv_halfZ]         = fabs(halez);
}

Acts::CylinderVolumeBounds::CylinderVolumeBounds(double rinner, double router, double haphi, double halez) :
 VolumeBounds(),
 m_boundValues(4,0.)
{
    m_boundValues[bv_innerRadius]   = fabs(rinner);
    m_boundValues[bv_outerRadius]   = fabs(router);
    m_boundValues[bv_halfPhiSector] = fabs(haphi);
    m_boundValues[bv_halfZ]         = fabs(halez);
 }

Acts::CylinderVolumeBounds::CylinderVolumeBounds(const Acts::CylinderVolumeBounds& cylbo) :
 VolumeBounds(),
 m_boundValues(cylbo.m_boundValues)
{}

Acts::CylinderVolumeBounds::~CylinderVolumeBounds()
{}

Acts::CylinderVolumeBounds& Acts::CylinderVolumeBounds::operator=(const Acts::CylinderVolumeBounds& cylbo)
{
  if (this!=&cylbo)
      m_boundValues    = cylbo.m_boundValues;
  return *this;
}

const std::vector<const Acts::Surface*>* Acts::CylinderVolumeBounds::decomposeToSurfaces(std::shared_ptr<Acts::Transform3D> transformPtr) const
{
    std::vector<const Acts::Surface*>* retsf = new std::vector<const Acts::Surface*>;
    // memory optimisation --- reserve the maximum
    retsf->reserve(6);

    // set the transform
    Acts::Transform3D transform = ( transformPtr == nullptr) ? Acts::Transform3D::Identity() : (*transformPtr.get());
    Acts::Transform3D* tTransform = nullptr;
    Acts::RotationMatrix3D discRot(transform.rotation());
    Acts::Vector3D cylCenter(transform.translation());

    // bottom Disc (negative z)
    Acts::RotationMatrix3D bottomDiscRot;
    bottomDiscRot.col(0) = discRot.col(1);
    bottomDiscRot.col(1) = discRot.col(0);
    bottomDiscRot.col(2) = -discRot.col(2);

    std::shared_ptr<const Acts::RadialBounds> dBounds(discBounds());
    tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(M_PI, Acts::Vector3D(1.,0.,0.))*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
    retsf->push_back(new Acts::DiscSurface(std::shared_ptr<Acts::Transform3D>(tTransform),dBounds));
    // top Disc (positive z)
    tTransform = new Acts::Transform3D(discRot*Acts::Translation3D(cylCenter + halflengthZ()*discRot.col(2)));
    retsf->push_back(new Acts::DiscSurface(std::shared_ptr<Acts::Transform3D>(tTransform),dBounds));

    // outer Cylinder - shares the transform
    retsf->push_back(new Acts::CylinderSurface(transformPtr, outerCylinderBounds()));

    // innermost Cylinder
    if (innerRadius() > s_numericalStable )
      retsf->push_back(new Acts::CylinderSurface(transformPtr, innerCylinderBounds()));

    // the cylinder is sectoral
    if ( fabs(halfPhiSector() - M_PI) >s_numericalStable) {
      std::shared_ptr<const Acts::PlanarBounds> sp12Bounds(sectorPlaneBounds());
      // sectorPlane 1 (negative phi)
      Acts::Transform3D* sp1Transform = new Acts::Transform3D(transform*Acts::AngleAxis3D(-halfPhiSector(), Acts::Vector3D(0.,0.,1.))
                                           *Acts::Translation3D(Acts::Vector3D(mediumRadius(),0.,0.))*Acts::AngleAxis3D(M_PI/2,Acts::Vector3D(1.,0.,0.)));
      retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(sp1Transform),sp12Bounds));
      // sectorPlane 2 (positive phi)
      Acts::Transform3D* sp2Transform = new Acts::Transform3D(transform*Acts::AngleAxis3D(halfPhiSector(), Acts::Vector3D(0.,0.,1.))
                                           *Acts::Translation3D(Acts::Vector3D(mediumRadius(),0.,0.))*Acts::AngleAxis3D(-M_PI/2, Acts::Vector3D(1.,0.,0.)));
      retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(sp2Transform),sp12Bounds));
    }
    return retsf;
}


Acts::CylinderBounds* Acts::CylinderVolumeBounds::innerCylinderBounds() const
{
    return new Acts::CylinderBounds(m_boundValues[bv_innerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}

Acts::CylinderBounds* Acts::CylinderVolumeBounds::outerCylinderBounds() const
{
    return new Acts::CylinderBounds(m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}

Acts::RadialBounds* Acts::CylinderVolumeBounds::discBounds() const
{
    return new Acts::RadialBounds(m_boundValues[bv_innerRadius], m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector]);
}

Acts::RectangleBounds* Acts::CylinderVolumeBounds::sectorPlaneBounds() const
{
    return new Acts::RectangleBounds(0.5*(m_boundValues[bv_outerRadius]-m_boundValues[bv_innerRadius]), m_boundValues[bv_halfZ]);
}

std::ostream& Acts::CylinderVolumeBounds::dump( std::ostream& sl ) const
{
    return dumpT<std::ostream>(sl);
}
