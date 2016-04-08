///////////////////////////////////////////////////////////////////
// BevelledCylinderVolumeBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry  module
#include "Volumes/BevelledCylinderVolumeBounds.h"
#include "Surfaces/CylinderSurface.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/DiscSurface.h"
#include "Surfaces/SubtractedDiscSurface.h"
#include "Surfaces/EllipseBounds.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Surfaces/RectangleBounds.h"
#include "Surfaces/RadialBounds.h"
#include "Volumes/CuboidVolumeBounds.h"
#include "Volumes/CombinedVolumeBounds.h"
#include "Volumes/VolumeExcluder.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SystemOfUnits.h"
// STD/STL
#include <iostream>
#include <math.h>

double Acts::BevelledCylinderVolumeBounds::s_numericalStable = 10e-2 * Gaudi::Units::millimeter;

Acts::BevelledCylinderVolumeBounds::BevelledCylinderVolumeBounds() :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{}

Acts::BevelledCylinderVolumeBounds::BevelledCylinderVolumeBounds(double rinner, double router, double haphi, double halez, int type) :
 VolumeBounds(),
 m_boundValues(bv_length,0.),
 m_type(type),
 m_subtractedVolume(nullptr) 
 
{
    m_boundValues[bv_innerRadius]   = fabs(rinner);
    m_boundValues[bv_outerRadius]   = fabs(router);
    m_boundValues[bv_halfPhiSector] = fabs(haphi);
    m_boundValues[bv_halfZ]         = fabs(halez);
    m_boundValues[bv_thetaMinus]    = 0.;
    m_boundValues[bv_thetaPlus]     = 0.;   
}

Acts::BevelledCylinderVolumeBounds::BevelledCylinderVolumeBounds(const Acts::BevelledCylinderVolumeBounds& cylbo) :
 VolumeBounds(),
 m_boundValues(cylbo.m_boundValues),
 m_type(cylbo.m_type), 
 m_subtractedVolume(nullptr)
{}

Acts::BevelledCylinderVolumeBounds::~BevelledCylinderVolumeBounds()
{
  delete m_subtractedVolume;
}

Acts::BevelledCylinderVolumeBounds& 
   Acts::BevelledCylinderVolumeBounds::operator=(
     const Acts::BevelledCylinderVolumeBounds& cylbo)
{
  if (this!=&cylbo){
      m_boundValues      = cylbo.m_boundValues;
      m_type             = cylbo.m_type;
      m_subtractedVolume = 0;     
  }
  return *this;
}


const std::vector<const Acts::Surface*>* Acts::BevelledCylinderVolumeBounds::decomposeToSurfaces(std::shared_ptr<Acts::Transform3D> transformPtr) const
{
    std::vector<const Acts::Surface*>* retsf = new std::vector<const Acts::Surface*>;
    
    // memory optimisation (reserve a save number of 20)
    retsf->reserve(6);
    
    // the transform
    Acts::Transform3D transform = ( transformPtr == nullptr) ? Acts::Transform3D::Identity() : (*transformPtr.get());
    Acts::RotationMatrix3D discRot(transform.rotation());
    Acts::Vector3D  cylCenter(transform.translation());

    if (m_type > -1 && !m_subtractedVolume) m_subtractedVolume = subtractedVolume();  
   
    // bottom Ellipse/Disc (negative z)
    Acts::Transform3D* tTransform = nullptr;
                                        
    if (m_type<0) {
        tTransform = new Acts::Transform3D((discRot*Acts::AngleAxis3D(-m_boundValues[bv_thetaMinus]+M_PI, Acts::Vector3D(0.,1.,0.))) *
                                               Acts::Translation3D(cylCenter - (halflengthZ()-m_boundValues[bv_outerRadius]*tan(m_boundValues[bv_thetaMinus]))*discRot.col(2)));
        retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),bottomEllipseBounds()));
    } else {
      if (m_subtractedVolume) {
	      // create the subracted volume
          Acts::Volume* subtrVol = new Acts::Volume(*m_subtractedVolume); 
          tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(M_PI, Acts::Vector3D(1.,0.,0.))*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
	      Acts::DiscSurface bottomDisc(std::shared_ptr<Acts::Transform3D>(tTransform),discBounds());
	      retsf->push_back(new Acts::SubtractedDiscSurface(bottomDisc,new Acts::VolumeExcluder(subtrVol),false));
      } else {
          tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(M_PI, Acts::Vector3D(1.,0.,0.))*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
          retsf->push_back(new Acts::DiscSurface(std::shared_ptr<Acts::Transform3D>(tTransform),discBounds()));   
       }
    }

    // top Ellipse/Disc (positive z)
    if (m_type<0) {
        tTransform = new Acts::Transform3D(discRot*Acts::AngleAxis3D(m_boundValues[bv_thetaPlus], Acts::Vector3D(0.,1.,0.))*
                     Acts::Translation3D(cylCenter + (halflengthZ()-m_boundValues[bv_outerRadius]*tan(m_boundValues[bv_thetaPlus]))*discRot.col(2))); 
        retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),topEllipseBounds())); 
    } else {
      if (m_subtractedVolume) {
	      Acts::Volume* subtrVol = new Acts::Volume(*m_subtractedVolume); 
          tTransform = new Acts::Transform3D(transform*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
	      Acts::DiscSurface topDisc(std::shared_ptr<Acts::Transform3D>(tTransform), discBounds());
	      retsf->push_back(new Acts::SubtractedDiscSurface(topDisc,new Acts::VolumeExcluder(subtrVol),false));
      } else {
          tTransform = new Acts::Transform3D(transform*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
          retsf->push_back(new Acts::DiscSurface(std::shared_ptr<Acts::Transform3D>(tTransform), discBounds()));
      }
    }

    // outer BevelledCylinder/Plane
    if (m_type<0) 
        retsf->push_back(new Acts::CylinderSurface(transformPtr, outerBevelledCylinderBounds()));
    else if (m_type<2) 
        retsf->push_back(new Acts::CylinderSurface(transformPtr,outerCylinderBounds()));
    else { 
        tTransform = new Acts::Transform3D(transform *Acts::Translation3D(Acts::Vector3D(this->outerRadius(),0.,0.))
                                        * Acts::AngleAxis3D(-90 *Gaudi::Units::deg, Acts::Vector3D(0.,0.,1.))*Acts::AngleAxis3D(-90.*Gaudi::Units::deg, Acts::Vector3D(1.,0.,0.)));
        retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),outerBevelledPlaneBounds()));
    }

    // inner BevelledCylinder/Plane
    if (innerRadius() > s_numericalStable ){
      if (m_type<1) 
          retsf->push_back(new Acts::CylinderSurface(transformPtr, innerBevelledCylinderBounds()));
      else if (m_type==2) 
          retsf->push_back(new Acts::CylinderSurface(transformPtr,innerCylinderBounds()));
      else {
          tTransform = new Acts::Transform3D(transform *Acts::Translation3D(Acts::Vector3D(this->innerRadius(),0.,0.))
						  *Acts::AngleAxis3D(+90*Gaudi::Units::deg, Acts::Vector3D(0.,0.,1.))*Acts::AngleAxis3D(-90.*Gaudi::Units::deg, Acts::Vector3D(1.,0.,0.)));
          retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),innerBevelledPlaneBounds()));
    
        }
    }         
    
    // phi-sectoral
    if ( fabs(halfPhiSector() - M_PI) >s_numericalStable) {
      if ( m_type<0 ) {
	    // sectorPlane 1 (negative phi)     
        std::shared_ptr<const Acts::PlanarBounds> trdBounds(sectorTrdBounds());      
        tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(-halfPhiSector(), Acts::Vector3D(0.,0.,1.))
	    							  *Acts::Translation3D(Acts::Vector3D(mediumRadius(),0.,0.))*Acts::AngleAxis3D(M_PI/2, Acts::Vector3D(1.,0.,0.)));      
	    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),trdBounds));        
	    // sectorPlane 2 (positive phi)    
        tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(halfPhiSector(), Acts::Vector3D(0.,0.,1.))
	    							  *Acts::Translation3D(Acts::Vector3D(mediumRadius(),0.,0.))*Acts::AngleAxis3D(-M_PI/2, Acts::Vector3D(1.,0.,0.)));
	    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),trdBounds));    
      } else {
	    // sectorPlane 1 (negative phi)    
	    double ri = innerRadius(); double ro = outerRadius();
	    if (m_type==1 || m_type==3) ri *= 1./cos(halfPhiSector()); 
	    if (m_type>1 ) ro *= 1./cos(halfPhiSector()); 
        std::shared_ptr<const Acts::PlanarBounds> sBounds(sectorPlaneBounds());
        tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(-halfPhiSector(), Acts::Vector3D(0.,0.,1.))
	    							  *Acts::Translation3D(Acts::Vector3D(0.5*(ri+ro),0.,0.))*Acts::AngleAxis3D(M_PI/2, Acts::Vector3D(1.,0.,0.)));
	    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),sBounds));        
	    // sectorPlane 2 (positive phi)
        tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(halfPhiSector(), Acts::Vector3D(0.,0.,1.))
	    							  *Acts::Translation3D(Acts::Vector3D(0.5*(ri+ro),0.,0.))*Acts::AngleAxis3D(-M_PI/2, Acts::Vector3D(1.,0.,0.)));
	    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),sBounds));
      }    
    }
    return retsf;
}
    
Acts::CylinderBounds* Acts::BevelledCylinderVolumeBounds::innerBevelledCylinderBounds() const
{

  return new Acts::CylinderBounds(m_boundValues[bv_innerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}
    
Acts::CylinderBounds* Acts::BevelledCylinderVolumeBounds::outerBevelledCylinderBounds() const
{
  return new Acts::CylinderBounds(m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}

Acts::RectangleBounds* Acts::BevelledCylinderVolumeBounds::outerBevelledPlaneBounds() const
{
    return new Acts::RectangleBounds(m_boundValues[bv_outerRadius]*tan(m_boundValues[bv_halfPhiSector]), m_boundValues[bv_halfZ]);
}

Acts::RectangleBounds* Acts::BevelledCylinderVolumeBounds::innerBevelledPlaneBounds() const
{
    return new Acts::RectangleBounds(m_boundValues[bv_innerRadius]*tan(m_boundValues[bv_halfPhiSector]), m_boundValues[bv_halfZ]);
}
    
Acts::EllipseBounds* Acts::BevelledCylinderVolumeBounds::bottomEllipseBounds() const
{
//    return new Acts::EllipseBounds();
    return new Acts::EllipseBounds(m_boundValues[bv_innerRadius]/cos(m_boundValues[bv_thetaMinus]),m_boundValues[bv_innerRadius], m_boundValues[bv_outerRadius]/cos(m_boundValues[bv_thetaMinus]), m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector]);
}

Acts::EllipseBounds* Acts::BevelledCylinderVolumeBounds::topEllipseBounds() const
{
//    return new Acts::EllipseBounds();
    return new Acts::EllipseBounds(m_boundValues[bv_innerRadius]/cos(m_boundValues[bv_thetaPlus]),m_boundValues[bv_innerRadius], m_boundValues[bv_outerRadius]/cos(m_boundValues[bv_thetaPlus]), m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector]);
    
}

Acts::CylinderBounds* Acts::BevelledCylinderVolumeBounds::innerCylinderBounds() const
{
    return new Acts::CylinderBounds(m_boundValues[bv_innerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}
    
Acts::CylinderBounds* Acts::BevelledCylinderVolumeBounds::outerCylinderBounds() const
{
    return new Acts::CylinderBounds(m_boundValues[bv_outerRadius], m_boundValues[bv_halfPhiSector], m_boundValues[bv_halfZ]);
}
    
Acts::RadialBounds* Acts::BevelledCylinderVolumeBounds::discBounds() const
{
  // adjust radius to make sure all surface covered
  double outerRadius = ( m_type>1 ) ? m_boundValues[bv_outerRadius]/cos(m_boundValues[bv_halfPhiSector]) : m_boundValues[bv_outerRadius]; 
  return new Acts::RadialBounds(m_boundValues[bv_innerRadius], outerRadius, m_boundValues[bv_halfPhiSector]);
}

Acts::TrapezoidBounds* Acts::BevelledCylinderVolumeBounds::sectorTrdBounds() const
{
  return new Acts::TrapezoidBounds(0.5*(outerRadius()-innerRadius()), m_boundValues[bv_halfZ], m_boundValues[bv_thetaMinus], m_boundValues[bv_thetaPlus]);
}

Acts::RectangleBounds* Acts::BevelledCylinderVolumeBounds::sectorPlaneBounds() const
{
  double ri = innerRadius(); double ro = outerRadius();
  if (m_type==1 || m_type==3) ri *= 1./cos(halfPhiSector()); 
  if (m_type>1 ) ro *= 1./cos(halfPhiSector()); 
  return new Acts::RectangleBounds(0.5*(ro-ri), m_boundValues[bv_halfZ]);
}

Acts::Volume* Acts::BevelledCylinderVolumeBounds::subtractedVolume() const
{
  if (m_type<1) return 0; 

  double tp = tan(m_boundValues[bv_halfPhiSector]);
  Acts::Volume* volIn = 0; 
  Acts::Volume* volOut = 0; 
  if ( m_type==1 || m_type==3 ) {   // cut inner cylinder  
    volIn = new Acts::Volume(0,new Acts::CuboidVolumeBounds(m_boundValues[bv_innerRadius],m_boundValues[bv_innerRadius]*tp+0.1,m_boundValues[bv_halfZ] + 0.1) );
  }
  if ( m_type>1 ) { 
    double hz = m_boundValues[bv_outerRadius]*(1./cos(m_boundValues[bv_halfPhiSector])-1.); 
    volOut = new Acts::Volume(new Acts::Transform3D( Acts::Translation3D(Acts::Vector3D(m_boundValues[bv_outerRadius]+hz,0.,0.))),
			     new Acts::CuboidVolumeBounds(hz,m_boundValues[bv_outerRadius]*tp+0.1,m_boundValues[bv_halfZ] +0.1) );
  }
 
  if (!volIn)  m_subtractedVolume = volOut;
  else if (!volOut) m_subtractedVolume = volIn;
  else m_subtractedVolume = new Acts::Volume(0, new Acts::CombinedVolumeBounds(volIn,volOut,false));

  return m_subtractedVolume; 
}

// ostream operator overload
MsgStream& Acts::BevelledCylinderVolumeBounds::dump( MsgStream& sl ) const
{
    return dumpT<MsgStream>(sl);
}

std::ostream& Acts::BevelledCylinderVolumeBounds::dump( std::ostream& sl ) const 
{
    return dumpT<std::ostream>(sl);
}
