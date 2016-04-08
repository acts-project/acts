///////////////////////////////////////////////////////////////////
// DiscSurface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Trk
#include "Surfaces/DiscSurface.h"
#include "Surfaces/RadialBounds.h"
#include "Surfaces/DiscTrapezoidalBounds.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
// STD/STL
#include <iostream>
#include <iomanip>
//Amg
#include "Algebra/AlgebraDefinitions.h"

Acts::NoBounds Acts::DiscSurface::s_boundless;

// default constructor
Acts::DiscSurface::DiscSurface() :
  Acts::Surface(),
  m_bounds(nullptr)
{}

// copy constructor
Acts::DiscSurface::DiscSurface(const DiscSurface& dsf) :
  Acts::Surface(dsf),
  m_bounds(dsf.m_bounds)
{}

// copy constructor with shift
Acts::DiscSurface::DiscSurface(const DiscSurface& dsf, const Acts::Transform3D& transf ) :
  Acts::Surface(dsf, transf),
  m_bounds(dsf.m_bounds)
{}

// construct a disc with full phi coverage
Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, double rmin, double rmax) :
  Acts::Surface(htrans),
  m_bounds(new Acts::RadialBounds(rmin, rmax))
{}

// construct a disc with given phi coverage
Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, double rmin, double rmax, double hphisec) :
  Acts::Surface(htrans),
  m_bounds(new Acts::RadialBounds(rmin, rmax, hphisec))
{}

Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, double minhalfx, double maxhalfx, double maxR, double minR, double avephi, double stereo):
  Acts::Surface(htrans),
  m_bounds(new Acts::DiscTrapezoidalBounds(minhalfx, maxhalfx, maxR, minR, avephi, stereo)),
  m_referencePoint(nullptr)
{}

// construct a disc with given bounds
Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, const Acts::RadialBounds* dbounds) :
  Acts::Surface(htrans),
  m_bounds(dbounds),
  m_referencePoint(nullptr)
{}
        
// construct a disc with given bounds
Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, const Acts::DiscTrapezoidalBounds* dbounds) :
  Acts::Surface(htrans),
  m_bounds(dbounds)
{}
        
// construct a disc with given bounds
Acts::DiscSurface::DiscSurface(std::shared_ptr<Acts:: Transform3D> htrans, std::shared_ptr<const Acts::DiscBounds> dbounds) :
  Acts::Surface(htrans),
  m_bounds(dbounds)
{}

// construct a disc from a transform, bounds is not set.
Acts::DiscSurface::DiscSurface(std::unique_ptr<Acts::Transform3D> htrans) :
  Acts::Surface(std::shared_ptr<Acts::Transform3D>(std::move(htrans))),
  m_bounds(nullptr)
{}

// construct form DetectorElementBase
Acts::DiscSurface::DiscSurface(const Acts::DetectorElementBase& detelement, const Identifier& identifier) :
  Acts::Surface(detelement, identifier),
  m_bounds(nullptr)
{}

// destructor (will call destructor from base class which deletes objects)
Acts::DiscSurface::~DiscSurface()
{}

Acts::DiscSurface& Acts::DiscSurface::operator=(const DiscSurface& dsf)
{
  if (this!=&dsf){
    Acts::Surface::operator=(dsf);
    m_bounds =  dsf.m_bounds;
  }
  return *this;
}

bool Acts::DiscSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::DiscSurface* dsf = dynamic_cast<const Acts::DiscSurface*>(&sf);
  if (!dsf) return false;
  if (this==dsf) return true;
  bool transfEqual(transform().isApprox(dsf->transform(),10e-8));
  bool centerEqual = (transfEqual) ? (center() == dsf->center()) : false;
  bool boundsEqual = (centerEqual) ? (bounds() == dsf->bounds()) : false;
  return boundsEqual;
}

void Acts::DiscSurface::localToGlobal(const Acts::Vector2D& locpos, const Acts::Vector3D&, Acts::Vector3D& glopos) const
{ 
  // create the position in the local 3d frame
  Acts::Vector3D loc3Dframe(locpos[Acts::eLOC_R]*cos(locpos[Acts::eLOC_PHI]),
			   locpos[Acts::eLOC_R]*sin(locpos[Acts::eLOC_PHI]),
			   0.);
  // transport it to the globalframe
  glopos = transform()*loc3Dframe;
}

/** local<->global transformation in case of polar local coordinates */
bool Acts::DiscSurface::globalToLocal(const Acts::Vector3D& glopos, const Acts::Vector3D&, Acts::Vector2D& locpos) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse())*glopos;
  locpos = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.phi() );
  return ( ( fabs(loc3Dframe.z()) > s_onSurfaceTolerance ) ? false : true );
}

const Acts::Vector2D Acts::DiscSurface::localPolarToLocalCartesian(const Acts::Vector2D& locpol) const
{
  const Acts::DiscTrapezoidalBounds* dtbo = dynamic_cast<const Acts::DiscTrapezoidalBounds*>(&(bounds()));
  if (dtbo) {
    double rMedium = dtbo->rCenter();
    double phi     = dtbo->averagePhi();

    Acts::Vector2D polarCenter(rMedium, phi);
    Acts::Vector2D cartCenter = localPolarToCartesian(polarCenter);
    Acts::Vector2D cartPos = localPolarToCartesian(locpol);
    Acts::Vector2D Pos = cartPos - cartCenter;
    
    Acts::Vector2D locPos(Pos[Acts::eLOC_X]*sin(phi) - Pos[Acts::eLOC_Y]*cos(phi),
			 Pos[Acts::eLOC_Y]*sin(phi) + Pos[Acts::eLOC_X]*cos(phi)); 
    
    return Acts::Vector2D(locPos[Acts::eLOC_X], locPos[Acts::eLOC_Y]); 
  }
  return Acts::Vector2D(locpol[Acts::eLOC_R]*cos(locpol[Acts::eLOC_PHI]),locpol[Acts::eLOC_R]*sin(locpol[Acts::eLOC_PHI])); 
}

/** local<->global transformation in case of polar local coordinates */
const Acts::Vector3D Acts::DiscSurface::localCartesianToGlobal(const Acts::Vector2D& locpos) const 
{
  Acts::Vector3D loc3Dframe(locpos[Acts::eLOC_X], locpos[Acts::eLOC_Y], 0.);
  return Acts::Vector3D(transform()*loc3Dframe);
}

/** local<->global transformation in case of polar local coordinates */
const Acts::Vector2D Acts::DiscSurface::globalToLocalCartesian(const Acts::Vector3D& glopos, double) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse())*glopos;
  return Acts::Vector2D(loc3Dframe.x(), loc3Dframe.y());
}

bool Acts::DiscSurface::isOnSurface(const Acts::Vector3D& glopo, const BoundaryCheck& bchk) const
{
  Acts::Vector3D loc3Dframe = (transform().inverse())*glopo;
  if ( fabs(loc3Dframe.z()) > (s_onSurfaceTolerance) ) return false;
  return (bchk ? bounds().inside(Acts::Vector2D(loc3Dframe.perp(),loc3Dframe.phi()),bchk) : true);
}

