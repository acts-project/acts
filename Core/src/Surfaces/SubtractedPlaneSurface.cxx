///////////////////////////////////////////////////////////////////
// SubtractedPlaneSurface.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/SubtractedPlaneSurface.h"

// default constructor
Acts::SubtractedPlaneSurface::SubtractedPlaneSurface() :
  Acts::PlaneSurface(),
  m_subtrVol(),
  m_shared(true)
{}

// copy constructor
Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(const SubtractedPlaneSurface& psf) :
  Acts::PlaneSurface(psf),
  m_subtrVol(psf.m_subtrVol),
  m_shared(psf.m_shared)
{}

// copy constructor with shift
Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(const SubtractedPlaneSurface& psf, const Acts::Transform3D& transf) :
  Acts::PlaneSurface(psf, transf),
  m_subtrVol(psf.m_subtrVol),
  m_shared(psf.m_shared)
{}

// constructor
Acts::SubtractedPlaneSurface::SubtractedPlaneSurface(const Acts::PlaneSurface& ps, AreaExcluder* vol, bool shared) :
  Acts::PlaneSurface(ps),
  m_subtrVol(vol),
  m_shared(shared)
{}

// destructor (will call destructor from base class which deletes objects)
Acts::SubtractedPlaneSurface::~SubtractedPlaneSurface()
{}

Acts::SubtractedPlaneSurface& Acts::SubtractedPlaneSurface::operator=(const Acts::SubtractedPlaneSurface& psf){
  
  if (this!=&psf){
    Acts::PlaneSurface::operator=(psf);
    m_subtrVol = psf.m_subtrVol;
    m_shared = psf.m_shared;
  }
  return *this;

} 

bool Acts::SubtractedPlaneSurface::operator==(const Acts::Surface& sf) const
{
  // first check the type not to compare apples with oranges
  const Acts::SubtractedPlaneSurface* spsf = dynamic_cast<const Acts::SubtractedPlaneSurface*>(&sf);
  if (!spsf) return false;
    bool surfaceEqual = Acts::PlaneSurface::operator==(sf);
    bool sharedEqual = (surfaceEqual) ? (shared() == spsf->shared()) : false; 
    return sharedEqual;
}


