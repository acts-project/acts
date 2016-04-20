///////////////////////////////////////////////////////////////////
// SlidingDiscSurface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Surfaces/SlidingDiscSurface.h"
#include "ACTS/Surfaces/DiscBounds.h"
#include "ACTS/Utilities/PrecisionDefinition.h"

// default constructor
Acts::SlidingDiscSurface::SlidingDiscSurface() :
  Acts::DiscSurface(),
  m_depth(),
  m_etaBin(),
  m_align()
{}

// copy constructor
Acts::SlidingDiscSurface::SlidingDiscSurface(const SlidingDiscSurface& dsf) :
  Acts::DiscSurface(dsf),
  m_depth(new std::vector<float>(*(dsf.m_depth))),
  m_etaBin(dsf.m_etaBin->clone()),
  m_align(dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr)
{}

// copy constructor with shift
Acts::SlidingDiscSurface::SlidingDiscSurface(const SlidingDiscSurface& dsf, const Acts::Transform3D& transf ) :
  Acts::DiscSurface(dsf, transf),
  m_depth(new std::vector<float>(*(dsf.m_depth))),
  m_etaBin(dsf.m_etaBin->clone()),
  m_align(dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr)
{}

// constructor
Acts::SlidingDiscSurface::SlidingDiscSurface(const Acts::DiscSurface& dsf, 
                                            Acts::BinUtility* bu, 
					                        const std::vector<float>* offset,
                                            Acts::Transform3D* align) :
  Acts::DiscSurface(dsf),
  m_depth(offset),
  m_etaBin(bu),
  m_align(align)
{}

// destructor (will call destructor from base class which deletes objects)
Acts::SlidingDiscSurface::~SlidingDiscSurface()
{
  delete m_depth;
  delete m_etaBin;
  delete m_align;
}

Acts::SlidingDiscSurface& Acts::SlidingDiscSurface::operator=(const SlidingDiscSurface& dsf)
{
  if (this!=&dsf){
   Acts::DiscSurface::operator=(dsf);
   delete m_depth;
   m_depth = new std::vector<float>(*(dsf.m_depth));
   delete m_etaBin;
   m_etaBin =  dsf.m_etaBin->clone();
   delete m_align;
   m_align = dsf.m_align ? new Acts::Transform3D(*dsf.m_align) : nullptr;
  }
  return *this;
}

bool Acts::SlidingDiscSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::SlidingDiscSurface* dsf = dynamic_cast<const Acts::SlidingDiscSurface*>(&sf);
  if (!dsf) return false;
  if (this==dsf) return true;
  bool transfEqual(transform().isApprox(dsf->transform(),10e-8));
  bool centerEqual = (transfEqual) ? (center() == dsf->center()) : false;
  bool boundsEqual = (centerEqual) ? (bounds() == dsf->bounds()) : false;
  return boundsEqual;
}

void Acts::SlidingDiscSurface::localToGlobal(const Acts::Vector2D& locpos, const Acts::Vector3D&, Acts::Vector3D& glopos) const
{
  // create the position in the local 3d frame
  Acts::Vector3D loc3D0(locpos[Acts::eLOC_R]*cos(locpos[Acts::eLOC_PHI]),
		       locpos[Acts::eLOC_R]*sin(locpos[Acts::eLOC_PHI]),
		       0.);
  // correct for alignment, retrieve offset correction
  Acts::Transform3D t0 = m_align ? m_align->inverse()*transform() : transform();
  float offset =  m_depth ? (*m_depth)[m_etaBin->bin(t0*loc3D0)] : 0.;
  Acts::Vector3D loc3Dframe(locpos[Acts::eLOC_R]*cos(locpos[Acts::eLOC_PHI]),
			   locpos[Acts::eLOC_R]*sin(locpos[Acts::eLOC_PHI]),
			   offset);
  // transport it to the globalframe
  glopos = Acts::Surface::transform()*loc3Dframe;
}

/** local<->global transformation in case of polar local coordinates */
bool Acts::SlidingDiscSurface::globalToLocal(const Acts::Vector3D& glopos, const Acts::Vector3D&, Acts::Vector2D& locpos) const
{
  Acts::Vector3D loc3D0 = m_align ? m_align->inverse()*glopos : glopos;       // used to retrieve localEta bin
  Acts::Vector3D loc3Dframe(Acts::Surface::transform().inverse()*glopos);
  locpos = Acts::Vector2D(loc3Dframe.perp(), loc3Dframe.phi() );
  return ( ( fabs(loc3Dframe.z()-(*m_depth)[m_etaBin->bin(loc3D0)]) > s_onSurfaceTolerance ) ? false : true );
}

bool Acts::SlidingDiscSurface::isOnSurface(const Acts::Vector3D& glopo, const Acts::BoundaryCheck& bchk) const
{
  Acts::Vector3D loc3D0 = m_align ? m_align->inverse()*glopo : glopo;     // used to retrieve localEta bin
  Acts::Vector3D loc3Dframe = (transform().inverse())*glopo;
  float offset = (*m_depth)[m_etaBin->bin(loc3D0)];
  if ( fabs(loc3Dframe.z()-offset) > (s_onSurfaceTolerance) ) return false;
  return (bchk ? bounds().inside(Acts::Vector2D(loc3Dframe.perp(),loc3Dframe.phi()), bchk) : true);
}
