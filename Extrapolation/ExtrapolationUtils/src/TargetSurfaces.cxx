///////////////////////////////////////////////////////////////////
//   Implementation file for class Ats::TargetSurfaces
///////////////////////////////////////////////////////////////////
// ATS project
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
// Version 1.0 09/2015 sarka.todorova@cern.ch  
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#include "GaudiKernel/MsgStream.h"
#include "ExtrapolationUtils/TargetSurfaces.h"
#include "Surfaces/Surface.h"
#include "Detector/TrackingVolume.h"

Ats::ExtrapolationCode  Ats::TargetSurfaces::setOnInput(Ats::ExCellCharged eCell, const Ats::Surface* sf, 
							Ats::BoundaryCheck bcheck) const
{
  Ats::Vector3D gp  = eCell.leadParameters->position();	
  Ats::Vector3D mom = eCell.leadParameters->momentum().normalized();	
  Ats::PropDirection dir = eCell.propDirection;
  
  return setOnInput( gp, mom*dir, eCell.leadVolume, sf, bcheck);
}

  Ats::ExtrapolationCode  Ats::TargetSurfaces::setOnInput(Ats::Vector3D position, Ats::Vector3D direction,
							 const Ats::TrackingVolume* fVol, 
							  const Ats::Surface* sf, Ats::BoundaryCheck bcheck) const
{
  // clear cache
  m_baseSurfaces.clear();
  m_tempSurfaces.clear();
  m_ordered.clear();
  m_orderTrue = false;

  // destination surface(s) first
  if (sf) { 
    Ats::TargetSurface tt(sf,bcheck,Ats::SurfNavigType::Target,0,0,Ats::TVNavigType::Unknown);
    evaluateInputDistance(tt,position,direction,true);  
    // abort if no valid intersection
    if (!m_baseSurfaces.size() || (m_baseSurfaces.back().distanceAlongPath<0 && !sf->isOnSurface(position,bcheck,m_tolerance)))
      return Ats::ExtrapolationCode::FailureDestination;
  }

  if ( initFrameVolume(position,direction,fVol) ) return Ats::ExtrapolationCode::InProgress;
  else   return Ats::ExtrapolationCode::FailureLoop;   // failure navigation?

}

Ats::TargetSurfaceVector  Ats::TargetSurfaces::orderedIntersections(Ats::ExCellNeutral eCell, const Ats::Surface* sf,
							       Ats::BoundaryCheck bcheck) const
{
  Ats::Vector3D gp  = eCell.leadParameters->position();	
  Ats::Vector3D mom = eCell.leadParameters->momentum().normalized();	
  Ats::PropDirection dir = eCell.propDirection;
  
  return orderedIntersections( gp, mom*dir, eCell.leadVolume, sf, bcheck);
}

Ats::TargetSurfaceVector  Ats::TargetSurfaces::orderedIntersections(Ats::Vector3D position, Ats::Vector3D direction,
					  const Ats::TrackingVolume* fVol, const Ats::Surface* sf,Ats::BoundaryCheck bcheck  ) const
{
  Ats::TargetSurfaceVector empty;   

  // clear cache
  m_baseSurfaces.clear();
  m_tempSurfaces.clear();
  m_ordered.clear();
  m_orderTrue = true;

  // destination surface(s) first
  if (sf) { 
    Ats::TargetSurface tt(sf,bcheck,Ats::SurfNavigType::Target,0,0,Ats::TVNavigType::Unknown);
    evaluateInputDistance(tt,position,direction,true);  
    // abort if no valid intersection
    if (!m_baseSurfaces.size() || (m_baseSurfaces.back().distanceAlongPath<0 && !sf->isOnSurface(position,bcheck,m_tolerance)))
      return empty;
  }

  if ( initFrameVolume(position,direction,fVol) ) return orderIntersections();
  else return empty;   // failure navigation?

}

bool Ats::TargetSurfaces::initFrameVolume(Ats::Vector3D pos, Ats::Vector3D dir, const Ats::TrackingVolume* fVol) const 
{
  if (!fVol) return true;

  m_currentFrame = fVol;             //  
  m_currentDense = fVol;             // default 

  if (m_debugMode) std::cout <<"DEBUG:input frame volume:"<<fVol->volumeName()<<" at position z:"<< fVol->center().z()<< std::endl;

  m_tolerance = 0.001; 

  // save probe coordinates
  m_probePos = pos;
  m_probeDir = dir; 

  // clear cache : keep target surfaces
  std::vector<Ats::TargetSurface>::iterator is = m_baseSurfaces.begin();
  while ( is!=m_baseSurfaces.end() && (*is).sfType==Ats::SurfNavigType::Target) is++;
  if (is!=m_baseSurfaces.end()) m_baseSurfaces.erase(is,m_baseSurfaces.end());   
  m_tempSurfaces.clear();

  // static frame boundaries
  auto& bounds = fVol->boundarySurfaces();
  for (unsigned int ib=0; ib< bounds.size(); ib++ ){
    const Ats::Surface& surf = (bounds[ib].get())->surfaceRepresentation();
    Ats::TargetSurface bb(&surf,true,Ats::SurfNavigType::BoundaryFrame,ib,fVol,Ats::TVNavigType::Frame);
    evaluateInputDistance(bb,pos,dir,true);
    if (m_debugMode) std::cout<<"DEBUG:frame input:id:status:distance:"<<ib<<":"<<bb.status<<":"<<bb.distanceAlongPath<<std::endl;
  }

  // check early exit
  m_nextSf = -1;
  double dExit = 0.;

  // index running over frame boundary only
  int index = 0;

  is=m_baseSurfaces.begin();
  while ( is!=m_baseSurfaces.end() ) {
    if ((*is).sfType==Ats::SurfNavigType::BoundaryFrame ) {
      double dist = (*is).distanceAlongPath;        
      if ( dist < m_tolerance && dist > dExit ) {
        Ats::Vector3D probe = pos+dir*dist;
        if ( (*is).surf->isOnSurface(probe,true,m_tolerance,m_tolerance) ) {
	  const Ats::TrackingVolume* nextVolume = bounds[(*is).index].get()->attachedVolume(probe,dir,Ats::alongMomentum); 
          if (nextVolume!=fVol) {
	    dExit = dist;
	    m_nextSf = index; 
	  }
        }  
      }
      index++;
    }
    is++;
  } 

  //if (m_debugMode) {
  //  findNext();
  //  std::cout <<"resolving SL exit from volume:"<<fVol->volumeName()<<":"<< m_nextSf<< ","<<dExit<<std::endl;
  //}

  if (m_nextSf>-1) {

    m_distanceToNext = dExit;
    if (m_debugMode) std::cout <<"DEBUG:early exit detected at boundary index:"<< m_nextSf<<","<<dExit <<std::endl;
    return false;
  }

  // fill frame volume

  // find closest surface and distance
  findNext();

  if (m_debugMode) std::cout<<"DEBUG:volume exit resolved (SL estimate):"<<m_nextSf<<","<<m_distanceToNext<< std::endl;

  return true;

}

void Ats::TargetSurfaces::evaluateInputDistance(Ats::TargetSurface& tt, Ats::Vector3D pos, Ats::Vector3D mom, bool base) const
{
    Ats::DistanceSolution distSol = tt.surf->straightLineDistanceEstimate(pos,mom);

    double dist = distSol.first();
    Ats::Vector3D posi =  pos + dist*mom;
    // skip trivial solutions
    if (distSol.numberOfSolutions()>1 && dist<m_tolerance && distSol.second()>m_tolerance) {
      dist = distSol.second();
      posi =  pos + dist*mom;
      if (m_orderTrue && !tt.surf->isOnSurface(posi,tt.bcheck,m_tolerance,m_tolerance) ) return;
      double dAbs = distSol.currentDistance(true);
      tt.setDistance(dist,fabs(dAbs),distSol.signedDistance() && dAbs!=0. ?  dAbs/fabs(dAbs)  : 0.);
      tt.setPosition(posi);
      if (fabs(dAbs)<m_tolerance) tt.setStatus(-1);
      else tt.setStatus(1);
      save(tt,base);
      return;
    } 
    // save closest solution 
    if (!m_orderTrue || dist>m_tolerance ) {
      double dAbs = distSol.currentDistance(true);
      tt.setDistance(dist,fabs(dAbs),distSol.signedDistance() && dAbs!=0. ?  dAbs/fabs(dAbs)  : 0.);
      tt.setPosition(posi);
      if (fabs(dAbs) < m_tolerance) tt.setStatus(-1);
      else if (dist > m_tolerance) tt.setStatus(1);
      save(tt,base);
    } 

    // save multiple intersection for neutral transport
    if (distSol.numberOfSolutions()>1 && distSol.second()>m_tolerance && m_orderTrue) {
      dist = distSol.second();
      posi =  pos + dist*mom;
      if ( tt.surf->isOnSurface(posi,tt.bcheck,m_tolerance,m_tolerance) ) {
        double dAbs = distSol.currentDistance(true);
        tt.setDistance(dist,fabs(dAbs),distSol.signedDistance() && dAbs!=0. ?  dAbs/fabs(dAbs)  : 0.);
        tt.setPosition(posi);
        if (fabs(dAbs)<m_tolerance) tt.setStatus(-1);
        else tt.setStatus(1);
        save(tt,base);
      }
    }

    return;
}  

bool Ats::TargetSurfaces::updateDistance(int index, Ats::TargetSurface& tt, Ats::Vector3D pos, Ats::Vector3D dir) const
{
    double previousDistance = tt.distanceAlongPath;

    Ats::DistanceSolution distSol = tt.surf->straightLineDistanceEstimate(pos,dir);

    double dist = 1.e08;

    // catch step across surface
    if (distSol.numberOfSolutions()>0 ) {
      dist = distSol.first();
      if (distSol.numberOfSolutions()>1 &&
      	  fabs(distSol.first()+m_lastStep-previousDistance) >
      	  fabs(distSol.second()+m_lastStep-previousDistance) && distSol.second() <= 0. ) {
	dist = distSol.second();  
      }
    }

    // verify true intersection and flip direction if necessary 
    if (  previousDistance*dist < 0. &&  fabs(dist)>m_tolerance ) {
      // verify change of sign in signedDistance ( after eliminating situations where this is meaningless )
      if ( !distSol.signedDistance() || fabs(distSol.currentDistance(true))<m_tolerance || tt.distance<m_tolerance
      	   || tt.signAbsDist*distSol.currentDistance(true)<0) {   // true intersection
	//if ( !distSol.signedDistance() || tt.signAbsDist*distSol.currentDistance(true)<0) {   // true intersection
	if (index==m_nextSf) {
          if (m_debugMode) std::cout <<"DEBUG:flipping intersection:signed ? true dist:"<< distSol.signedDistance()<<":"
				     << tt.signAbsDist*distSol.currentDistance(true) <<std::endl; 
	  m_flipDirection = true;
          m_flipDistance = dist;
	}  else if ( tt.distance>m_tolerance && fabs(distSol.currentDistance(true))>m_tolerance ) {
	  // here we need to compare with distance from current closest
	  if ( index>m_nextSf ) {   // easy case, already calculated
	    if ( dist<(m_flipDistance-m_tolerance) && tt.status!=-1)  {
	      m_flipDirection = true;
              m_flipDistance = dist;
	      m_nextSf=index;
	    }
	  } else if (m_distanceToNext>0.) {             // set as nearest (if not iterating already), will be rewritten later
	    if (tt.status!=-1) {
	      m_flipDirection = true;
              m_flipDistance = dist;
	      m_nextSf = index;
	    }
	  }
	}
      }
    }

    // continue iteration if appropriate
    if ( index == m_nextSf && dist<0. && previousDistance<dist ) m_distanceToNext = dist;
	
    // save current distance to surface
    Ats::Vector3D posi =  pos + dist*dir;
    double dAbs = distSol.currentDistance(true);
    tt.setDistance(dist,fabs(dAbs),distSol.signedDistance() && dAbs!=0. ?  dAbs/fabs(dAbs)  : 0.);
    tt.setPosition(posi);
    if (tt.status==-1 && fabs(dAbs)>m_tolerance ) tt.status = dist>0 ? 1 : 0;

    return true;
}  

void Ats::TargetSurfaces::save(Ats::TargetSurface& tt, bool base) const {

  if (base) m_baseSurfaces.push_back(tt);
  else {
    if (!m_tempSurfaces.size() || tt.assocVol!=m_tempSurfaces.back()[0].assocVol) {
      Ats::TargetSurfaceVector local; local.push_back(tt); m_tempSurfaces.push_back(local); 
    } else   m_tempSurfaces.back().push_back(tt);
  }        
}

Ats::TargetSurfaceVector  Ats::TargetSurfaces::orderIntersections() const{

  Ats::TargetSurfaceVector tsOrd ;		

  // base surfaces
  if (!m_baseSurfaces.size()) return tsOrd;
      
  std::vector<unsigned int> sols(m_baseSurfaces.size());
  for (unsigned int i=0;i<m_baseSurfaces.size(); i++) { sols[i]=i; }

  unsigned int itest=1;
  while ( itest<sols.size() ) {
    if ( m_baseSurfaces[sols[itest]].distanceAlongPath < m_baseSurfaces[sols[itest-1]].distanceAlongPath ) {
      unsigned int iex = sols[itest-1];
      sols[itest-1]=sols[itest];
      sols[itest]=iex;
      itest=1;
    } else itest++; 
  }

  for (unsigned int i=0;i<m_baseSurfaces.size(); i++) tsOrd.push_back(m_baseSurfaces[sols[i]]); 
  
  return tsOrd;

}

void Ats::TargetSurfaces::findNext() const{

  m_nextSf = -1;
  m_distanceToNext = 1.e06;
  m_numAlongPath = 0;

  // index running over all selected surfaces
  int index = 0;

  std::vector<Ats::TargetSurface>::iterator is=m_baseSurfaces.begin();
  while ( is!=m_baseSurfaces.end() ) {
    if ((*is).status==-1 && (*is).distanceAlongPath>m_tolerance) {
      m_numAlongPath++;
      double dd = (*is).distanceAlongPath;
      if ( dd < m_distanceToNext ) {
	m_nextSf = index;
	m_distanceToNext = dd ;
      }
    }
    if ((*is).status!=-1 && (*is).distanceAlongPath>m_tolerance) {
      m_numAlongPath++;
      double dd = (*is).distanceAlongPath;
      if (dd>m_tolerance && dd<(*is).distance ) dd=(*is).distance;
      if ( dd < m_distanceToNext ) {
	m_nextSf = index;
	m_distanceToNext = dd ; 
      }
    }
    index++; is++;
  } 
  for (unsigned it=0; it<m_tempSurfaces.size(); it++) {
    is=m_tempSurfaces[it].begin();
    while ( is!=m_tempSurfaces[it].end()) {
      if ((*is).status!=-1 && (*is).distanceAlongPath>m_tolerance) {
	m_numAlongPath++; 
	double dd = (*is).distanceAlongPath;
	if (dd>m_tolerance && dd<(*is).distance ) dd=(*is).distance;
	if ( dd < m_distanceToNext ) {
	  m_nextSf = index;
	  m_distanceToNext = dd; 
	}
      }
      index++; is++;
    }
  } 
}


bool Ats::TargetSurfaces::checkDistance(Ats::Vector3D pos, Ats::Vector3D dir, double nextStep) const{

  m_lastStep = (pos-m_probePos).mag();

  if (!m_lastStep>0.) {
    if (m_debugMode) std::cout <<"DEBUG:check distance with 0 step:"<<"next,dist:"<<m_nextSf<<","<<m_distanceToNext<< std::endl;
    return true;
  } 

  // dont overwrite previous estimates before full loop finished
  // limit the number of reevaluations by using rough estimates

  int nextSfCandidate = -1;
  m_distanceToNext = 1.e08;
  m_numAlongPath = 0;
  m_flipDirection = false;
  m_flipDistance = 0.;

  // index running over all selected surfaces
  int index = 0;

  std::vector<Ats::TargetSurface>::iterator is=m_baseSurfaces.begin();
  while ( is!=m_baseSurfaces.end() ) {
    // reevaluate distance if : status = -1  ( leaving surface with possibility of re-entry,
    //                                         switch to 0 when abs.distance > tolerance )
    //                          abs.distance < 2*nextStep
    // otherwise just subtract distance stepped for rough estimate                                         

    if ( (*is).status==-1 || ( fabs((*is).distance)-m_lastStep ) < 2*fabs(nextStep) || m_absDist) {  
      if (m_lastStep>m_tolerance) (*is).setStatus(0);    
      updateDistance(index,(*is),pos,dir);    
    } else {
      (*is).fastDistanceUpdate(m_lastStep);
    }

    if (m_debugMode) std::cout <<"DEBUG:check distance:"<<index<<":"<<(*is).status<<","<<(*is).distanceAlongPath<<
	      ","<<(*is).distance<<","<<(*is).signAbsDist<<":flip dir,dist:"<<m_flipDirection<<","<<m_flipDistance<< std::endl;
    if ( (*is).status!=-1 || (*is).distanceAlongPath>m_tolerance ) {
      if ( (*is).distanceAlongPath>-m_tolerance) {
	m_numAlongPath++; 
	double dd = (*is).distanceAlongPath;
	if (dd>m_tolerance && dd<(*is).distance ) dd=(*is).distance;
	if ( dd < m_distanceToNext ) {
	  nextSfCandidate = index;
	  m_distanceToNext = dd;
	} 
      }
    }   

    index++; is++;
  } 

  m_absDist=false;

  // distanceAlongPath negative : switch to absolute distance
  if (nextSfCandidate<0 && m_distanceToNext>0.) {
    if (!m_currentFrame->inside(pos,0.001)) {
      std::cout <<"ERROR:frame volume left, aborting:"<<m_currentFrame->volumeName()<<std::endl;
      return false; 
    }
    index = 0;
    is=m_baseSurfaces.begin();
    while ( is!=m_baseSurfaces.end() ) {
      updateDistance(index,(*is),pos,dir);    
      if ( (*is).status!=-1 || (*is).distance>m_tolerance ) {
	if ( (*is).distance>-m_tolerance) {
	  m_numAlongPath++; 
	  double dd = (*is).distance;
	  if ( dd < m_distanceToNext ) {
	    nextSfCandidate = index;
	    m_distanceToNext = dd;
	  } 
	}
      }   
      index++; is++;
    } 
    if (m_debugMode) std::cout <<"DEBUG:closest frame estimate based on absolute distance:"<< nextSfCandidate<<":"<<m_distanceToNext
	      <<", inside frame volume?"<< m_currentFrame->inside(pos,0.)<<std::endl; 
    m_absDist = true;
  }

  for (unsigned it=0; it<m_tempSurfaces.size(); it++) {
    is=m_tempSurfaces[it].begin();

    while ( is!=m_tempSurfaces[it].end()) {

      if ( (*is).status==-1 || ( fabs((*is).distance)-m_lastStep ) < 2*fabs(nextStep) ) {     
	updateDistance(index,(*is),pos,dir);    
      } else {
	(*is).fastDistanceUpdate(m_lastStep);
      }

      if ( (*is).status!=-1  || (*is).distanceAlongPath>m_tolerance ) {
	if ( (*is).distanceAlongPath>-m_tolerance) {
	  m_numAlongPath++; 
	  double dd = (*is).distanceAlongPath;
	  if (dd>m_tolerance && dd<(*is).distance ) dd=(*is).distance;
	  if ( dd<m_distanceToNext ) {
	    nextSfCandidate = index;
	    m_distanceToNext = dd;
	  } 
	}
      }

      index++; is++;
    }
  } 

  // distanceToNext estimate reliable for distances below 2*|nextStep| only 
  if ( !m_flipDirection && nextSfCandidate != m_nextSf  ) m_nextSf = nextSfCandidate;  
  //if ( !m_flipDirection && nextSfCandidate != m_nextSf && m_distanceToNext<2*fabs(nextStep) ) m_nextSf = nextSfCandidate;  

  // flip direction
  if ( m_flipDirection ) m_distanceToNext = m_flipDistance;

  m_probePos = pos;
  m_probeDir = dir;

  if (m_debugMode) std::cout <<"DEBUG:check distance returns:next,dist:"<<m_nextSf<<","<<m_distanceToNext<< std::endl;

  return true;
}

void Ats::TargetSurfaces::fillSolutions(int hitSf, Ats::Vector3D gp, TargetSurfaceVector& solutions) const
{
  // index running over all selected surfaces
  int index = 0;

  if (m_debugMode) std::cout <<"fill solutions at R,z,phi:"<<gp.perp()<<","<<gp.z()<<","<<gp.phi()<< std::endl;

  std::vector<Ats::TargetSurface>::iterator is=m_baseSurfaces.begin();
  while ( is!=m_baseSurfaces.end() ) {

    if ( (*is).status!=-1 ) {

      if ( index==hitSf || fabs((*is).distanceAlongPath)<0.01 ) {   // onSurface && boundary check
        if (m_debugMode) std::cout <<"DEBUG: onSurface, without bcheck, twice tolerance:"<<
			   (*is).surf->isOnSurface(gp,(*is).bcheck ,m_tolerance,m_tolerance)<<","<<
			   (*is).surf->isOnSurface(gp,false,m_tolerance,m_tolerance)<<","<<
			   (*is).surf->isOnSurface(gp,(*is).bcheck,2*m_tolerance,2*m_tolerance)<<std::endl;
	if ( (*is).surf->isOnSurface(gp,(*is).bcheck ,m_tolerance,m_tolerance) ) solutions.push_back((*is));
        else if ( (*is).surf->isOnSurface(gp,false,m_tolerance,m_tolerance) ) {   // probably numerical precision problem 
	   solutions.push_back((*is));                                            // in boundary check, 
           if (m_debugMode && !m_currentFrame->inside(gp,m_tolerance) )           // possible navigation break
	     std::cout<<"DEBUG: frame boundary crossed outside volume "<< m_currentFrame->volumeName()<< std::endl;  
	} else if ( index==hitSf) (*is).status = -1;
      }
    }     

    index++; is++;
  } 

  for (unsigned it=0; it<m_tempSurfaces.size(); it++) {
    is=m_tempSurfaces[it].begin();

    while ( is!=m_tempSurfaces[it].end()) {

      if ( (*is).status!=-1 ) {

	if ( index==hitSf || (*is).distanceAlongPath<0.01 ) {  // onSurface && boundary check
	  if ( (*is).surf->isOnSurface(gp,(*is).bcheck ,m_tolerance,m_tolerance) )  solutions.push_back((*is));
	  else if (index==hitSf) (*is).status=-1.;
	}
      }      

      index++; is++;
    }
  } 

  return;
}
