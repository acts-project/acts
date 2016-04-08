///////////////////////////////////////////////////////////////////
// TrapezoidVolumeBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "Volumes/TrapezoidVolumeBounds.h"
#include "GeometryUtils/GeometryStatics.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Surfaces/RectangleBounds.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SystemOfUnits.h"
// STD/STL
#include <iostream>
#include <iomanip>
#include <math.h>

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds() :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex, double maxhalex, double haley, double halez) :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{

  m_boundValues[bv_minHalfX] = minhalex;
  m_boundValues[bv_maxHalfX] = maxhalex; 
  m_boundValues[bv_halfY]    = haley;
  m_boundValues[bv_halfZ]    = halez;
  m_boundValues[bv_alpha]    = atan((m_boundValues[bv_maxHalfX]-m_boundValues[bv_minHalfX])/2/m_boundValues[bv_halfY]) + 0.5*M_PI;
  m_boundValues[bv_beta]     = m_boundValues[bv_alpha]; 

}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(double minhalex, double haley, double halez, double alpha, double beta) :
 VolumeBounds(),
 m_boundValues(bv_length,0.)
{
  m_boundValues[bv_minHalfX] = minhalex;
  m_boundValues[bv_halfY]    = haley;
  m_boundValues[bv_halfZ]    = halez;
  m_boundValues[bv_alpha]    = alpha;
  m_boundValues[bv_beta]     = beta; 
  // now calculate the remaining max half X
  double gamma = (alpha > beta) ? (alpha - 0.5*M_PI) : (beta - 0.5*M_PI);
  m_boundValues[bv_maxHalfX] = minhalex + (2.*haley)*tan(gamma);
}

Acts::TrapezoidVolumeBounds::TrapezoidVolumeBounds(const Acts::TrapezoidVolumeBounds& trabo) :
 VolumeBounds(),
 m_boundValues(trabo.m_boundValues)
{}

Acts::TrapezoidVolumeBounds::~TrapezoidVolumeBounds()
{}

Acts::TrapezoidVolumeBounds& Acts::TrapezoidVolumeBounds::operator=(const Acts::TrapezoidVolumeBounds& trabo)
{
  if (this!=&trabo)
      m_boundValues = trabo.m_boundValues;
  return *this;
}

const std::vector<const Acts::Surface*>* Acts::TrapezoidVolumeBounds::decomposeToSurfaces(std::shared_ptr<Acts::Transform3D> transformPtr) const
{
    std::vector<const Acts::Surface*>* retsf = new std::vector<const Acts::Surface*>;  
    Acts::Transform3D transform   = ( transformPtr == nullptr) ? Acts::Transform3D::Identity() : (*transformPtr.get());
    Acts::Transform3D* tTransform = nullptr; 
    // face surfaces xy
    Acts::RotationMatrix3D trapezoidRotation(transform.rotation());
    Acts::Vector3D  trapezoidX(trapezoidRotation.col(0));
    Acts::Vector3D  trapezoidY(trapezoidRotation.col(1));
    Acts::Vector3D  trapezoidZ(trapezoidRotation.col(2));
    Acts::Vector3D  trapezoidCenter(transform.translation());
    
    //   (1) - at negative local z
    std::shared_ptr<const PlanarBounds> xytBounds(faceXYTrapezoidBounds());
    tTransform = new Acts::Transform3D(transform*Acts::AngleAxis3D(180*Gaudi::Units::deg, Acts::Vector3D(0.,1.,0.))
							                                     *Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
      
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),xytBounds));    
    //   (2) - at positive local z
    tTransform = new Acts::Transform3D(transform*Acts::Translation3D(Acts::Vector3D(0.,0.,halflengthZ())));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),xytBounds));    
    // face surfaces yz
    // transmute cyclical
    //   (3) - at point A, attached to alpha opening angle 
    Acts::Vector3D A(minHalflengthX(), halflengthY(), trapezoidCenter.z());
    Acts::RotationMatrix3D alphaZRotation = (s_idRotation * Acts::AngleAxis3D (alpha() - 0.5*M_PI, Acts::Vector3D(0.,0.,1.))).toRotationMatrix();
    //CLHEP::HepRotation  alphaRotation(alphaZRotation*trapezoidRotation);
    Acts::RotationMatrix3D  faceAlphaRotation;
    faceAlphaRotation.col(0) = alphaZRotation.col(1);
    faceAlphaRotation.col(1) = -alphaZRotation.col(2);
    faceAlphaRotation.col(2) = -alphaZRotation.col(0);
    RectangleBounds* faceAlphaBounds = faceAlphaRectangleBounds();
    // Acts::Vector3D   faceAlphaPosition(A+faceAlphaRotation.colX()*faceAlphaBounds->halflengthX());
    Acts::Vector3D   faceAlphaPosition0(-0.5*(minHalflengthX()+maxHalflengthX()),0.,0.);
    Acts::Vector3D faceAlphaPosition = transform*faceAlphaPosition0;
    tTransform = new Acts::Transform3D((trapezoidRotation*faceAlphaRotation) * Acts::Translation3D(faceAlphaPosition));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),faceAlphaBounds));
    //   (4) - at point B, attached to beta opening angle 
    Acts::Vector3D B(minHalflengthX(), -halflengthY(), trapezoidCenter.z());
    Acts::RotationMatrix3D betaZRotation = (s_idRotation * Acts::AngleAxis3D (-(beta() - 0.5*M_PI),Acts::Vector3D(0.,0.,1.))).toRotationMatrix();
    //CLHEP::HepRotation  betaRotation(betaZRotation*trapezoidRotation);
    Acts::RotationMatrix3D  faceBetaRotation;
    faceBetaRotation.col(0) = betaZRotation.col(1);
    faceBetaRotation.col(1) = betaZRotation.col(2);
    faceBetaRotation.col(2) = betaZRotation.col(0);
    RectangleBounds* faceBetaBounds = faceBetaRectangleBounds();
    // Acts::Vector3D   faceBetaPosition(B+faceBetaRotation.colX()*faceBetaBounds->halflengthX());
    Acts::Vector3D   faceBetaPosition0( 0.5*(minHalflengthX()+maxHalflengthX()),0.,0.);
    Acts::Vector3D faceBetaPosition = transform*faceBetaPosition0;
    tTransform = new Acts::Transform3D(trapezoidRotation*faceBetaRotation * Acts::Translation3D(faceBetaPosition));
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),faceBetaBounds));
    // face surfaces zx
    //   (5) - at negative local x 
    tTransform = new Acts::Transform3D( transform * Acts::AngleAxis3D(180.*Gaudi::Units::deg, Acts::Vector3D(1.,0.,0.))
								*Acts::Translation3D(Acts::Vector3D(0.,halflengthY(),0.))
								*Acts::AngleAxis3D(-90*Gaudi::Units::deg, Acts::Vector3D(0.,1.,0.))*Acts::AngleAxis3D(-90.*Gaudi::Units::deg, Acts::Vector3D(1.,0.,0.))  );  
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), faceZXRectangleBoundsBottom()));
    //   (6) - at positive local x
    tTransform = new Acts::Transform3D( transform *Acts::Translation3D(Acts::Vector3D(0.,halflengthY(), 0.))
								*Acts::AngleAxis3D(-90*Gaudi::Units::deg, Acts::Vector3D(0.,1.,0.))*Acts::AngleAxis3D(-90.*Gaudi::Units::deg, Acts::Vector3D(1.,0.,0.))  );
    retsf->push_back(new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform), faceZXRectangleBoundsTop())); 

    return retsf;
}
    
// faces in xy
Acts::TrapezoidBounds* Acts::TrapezoidVolumeBounds::faceXYTrapezoidBounds() const
{
    return new Acts::TrapezoidBounds(m_boundValues[bv_minHalfX], m_boundValues[bv_maxHalfX], m_boundValues[bv_halfY]);
}

Acts::RectangleBounds* Acts::TrapezoidVolumeBounds::faceAlphaRectangleBounds() const
{
    return new Acts::RectangleBounds(m_boundValues[bv_halfY]/cos(m_boundValues[bv_alpha]-0.5*M_PI), m_boundValues[bv_halfZ]);
}

Acts::RectangleBounds* Acts::TrapezoidVolumeBounds::faceBetaRectangleBounds() const
{
    return new Acts::RectangleBounds(m_boundValues[bv_halfY]/cos(m_boundValues[bv_beta]-0.5*M_PI), m_boundValues[bv_halfZ]);
}

Acts::RectangleBounds* Acts::TrapezoidVolumeBounds::faceZXRectangleBoundsBottom() const
{
    return new Acts::RectangleBounds(m_boundValues[bv_halfZ], m_boundValues[bv_minHalfX]);
}

Acts::RectangleBounds* Acts::TrapezoidVolumeBounds::faceZXRectangleBoundsTop() const
{
    //double delta = (m_boundValues[bv_alpha] < m_boundValues[bv_beta]) ? m_boundValues[bv_alpha] - M_PI/2. : m_boundValues[bv_beta] - M_PI/2.;     
    // return new Acts::RectangleBounds(m_boundValues[bv_halfZ], 0.5*(m_boundValues[bv_minHalfX]+m_boundValues[bv_minHalfX]+2.*m_boundValues[bv_halfY]/cos(delta)));
    return new Acts::RectangleBounds(m_boundValues[bv_halfZ], m_boundValues[bv_maxHalfX]);
}

bool Acts::TrapezoidVolumeBounds::inside(const Acts::Vector3D& pos, double tol) const
{
    if (fabs(pos.z()) > m_boundValues[bv_halfZ] + tol) return false;
    if (fabs(pos.y()) > m_boundValues[bv_halfY] + tol) return false;
    Acts::TrapezoidBounds* faceXYBounds = faceXYTrapezoidBounds();
    Acts::Vector2D locp(pos.x(), pos.y());
    bool inside(faceXYBounds->inside(locp, BoundaryCheck(true, true, tol, tol)) ); 
    delete faceXYBounds;
    return inside;
}

// ostream operator overload
MsgStream& Acts::TrapezoidVolumeBounds::dump( MsgStream& sl ) const
{
    return dumpT<MsgStream>(sl);
}

std::ostream& Acts::TrapezoidVolumeBounds::dump( std::ostream& sl ) const 
{
    return dumpT<std::ostream>(sl);
}
