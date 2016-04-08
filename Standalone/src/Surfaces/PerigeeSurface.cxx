/////////////////////////////////////////////////////////////////
// PerigeeSurface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Surfaces/PerigeeSurface.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
// STD/STL
#include <iostream>
#include <iomanip>

Acts::NoBounds Acts::PerigeeSurface::s_perigeeBounds;

Acts::PerigeeSurface::PerigeeSurface() :
    Surface(),
    m_lineDirection(nullptr)
{}

Acts::PerigeeSurface::PerigeeSurface(const Acts::Vector3D& gp):
    Surface(),
    m_lineDirection(nullptr)
{
    Surface::m_center               = new Acts::Vector3D(gp);
    Surface::m_transform            = std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D());
    (*(Surface::m_transform.get())) = Acts::Translation3D(gp.x(),gp.y(),gp.z());
}

Acts::PerigeeSurface::PerigeeSurface(std::shared_ptr<Acts::Transform3D> tTransform):
    Surface(tTransform),
    m_lineDirection(nullptr)
{}

Acts::PerigeeSurface::PerigeeSurface(std::unique_ptr<Acts::Transform3D> tTransform):
  Surface(std::shared_ptr<Acts::Transform3D>(std::move(tTransform))),
  m_lineDirection(nullptr)
{}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& pesf):
    Surface(pesf),
    m_lineDirection(nullptr)
{}

Acts::PerigeeSurface::PerigeeSurface(const PerigeeSurface& pesf, const Acts::Transform3D& shift):
    Surface(),
    m_lineDirection(nullptr)
{
    if (pesf.m_center)         Surface::m_center = new Acts::Vector3D(shift*(*pesf.m_center));
    if (pesf.m_transform)      Surface::m_transform = std::shared_ptr<Acts::Transform3D>(new Acts::Transform3D(shift*(*pesf.m_transform)));
}

Acts::PerigeeSurface::~PerigeeSurface()
{
    delete m_lineDirection;
}

// assignment operator
Acts::PerigeeSurface& Acts::PerigeeSurface::operator=(const Acts::PerigeeSurface& pesf)
{
  if (this!=&pesf){
       Acts::Surface::operator=(pesf);
       delete m_lineDirection;
       m_lineDirection = new Acts::Vector3D(lineDirection());

  }
  return *this;
}

bool Acts::PerigeeSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::PerigeeSurface* persf = dynamic_cast<const Acts::PerigeeSurface*>(&sf);
  if (!persf) return false;
  if (this==persf) return true;
  return (center() == persf->center());
}


// true local to global method/
void Acts::PerigeeSurface::localToGlobal(const Acts::Vector2D& locpos,
                         				const Acts::Vector3D& glomom,
                         				Acts::Vector3D& glopos) const
{
    // this is for a tilted perigee surface
    if (Surface::m_transform){
        // get the vector perpenticular to the momentum and the straw axis
        Acts::Vector3D radiusAxisGlobal(lineDirection().cross(glomom));
        Acts::Vector3D locZinGlobal = transform()*Acts::Vector3D(0.,0.,locpos[Acts::eLOC_Z]);
        // transform zPosition into global coordinates and add Acts::eLOC_R * radiusAxis
        glopos = Acts::Vector3D(locZinGlobal + locpos[Acts::eLOC_R]*radiusAxisGlobal.normalized());
    } else {
        double phi = glomom.phi();
        glopos[Acts::eX] = - locpos[Acts::eLOC_D0]*sin(phi);
        glopos[Acts::eY] =   locpos[Acts::eLOC_D0]*cos(phi);
        glopos[Acts::eZ] =   locpos[Acts::eLOC_Z0];
        glopos += center();
    }
}

// true global to local method
bool Acts::PerigeeSurface::globalToLocal(const Acts::Vector3D& glopos,
                                        const Acts::Vector3D& glomom,
                                        Acts::Vector2D& locpos) const
{
    Acts::Vector3D perPos = (transform().inverse())*glopos;
    double d0 = perPos.perp();
    double z0 = perPos.z();
    // decide the sign of d0
    d0 *= ((lineDirection().cross(glomom)).dot(perPos)<0.0) ? -1.0 : 1.0;
    locpos[Acts::eLOC_D0] = d0;
    locpos[Acts::eLOC_Z0] = z0;
    return true;
}

// return the measurement frame - this is the frame where the covariance is defined
const Acts::RotationMatrix3D Acts::PerigeeSurface::measurementFrame(const Acts::Vector3D&, const Acts::Vector3D& glomom) const
{
    Acts::RotationMatrix3D mFrame;
    // construct the measurement frame
    const Acts::Vector3D& measY = lineDirection();
    Acts::Vector3D measX(measY.cross(glomom).unit());
    Acts::Vector3D measDepth(measX.cross(measY));
    // assign the columnes
    mFrame.col(0) = measX;
    mFrame.col(1) = measY;
    mFrame.col(2) = measDepth;
    // return the rotation matrix
    return mFrame;
}

// overload of ostream operator
MsgStream& Acts::PerigeeSurface::dump( MsgStream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::PerigeeSurface:" << std::endl;
    sl << "     Center position  (x, y, z) = (" << center().x() << ", " << center().y() << ", " << center().z() << ")";
    sl << std::setprecision(-1);
    return sl;
}

// overload of ostream operator
std::ostream& Acts::PerigeeSurface::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::PerigeeSurface:" << std::endl;
    sl << "     Center position  (x, y, z) = (" << center().x() << ", " << center().y() << ", " << center().z() << ")";
    sl << std::setprecision(-1);
    return sl;
}
