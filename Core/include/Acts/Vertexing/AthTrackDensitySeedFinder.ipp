/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/*********************************************************************
          TrackDensitySeedFinder.cxx - Description in header file
*********************************************************************/

#include "TrkVertexSeedFinderTools/TrackDensitySeedFinder.h"

#include "TrkVertexSeedFinderUtils/SeedFinderParamDefs.h"

#include "TrkParameters/TrackParameters.h"
#include "TrkTrack/Track.h"
#include "TrkEventPrimitives/ParamDefs.h"

//Amg
#include "GeoPrimitives/GeoPrimitives.h"

#include <algorithm>

namespace Trk
{

  TrackDensitySeedFinder::TrackDensitySeedFinder(const std::string& t, const std::string& n, const IInterface*  p) : 
    base_class(t,n,p)
  {   
  }


  TrackDensitySeedFinder::~TrackDensitySeedFinder()
  {
  }


  StatusCode TrackDensitySeedFinder::initialize() 
  { 
    ATH_CHECK( m_densityEstimator.retrieve() );
    ATH_MSG_DEBUG("Initialize successful");
    return StatusCode::SUCCESS;
  }

  StatusCode TrackDensitySeedFinder::finalize() 
  {
    ATH_MSG_DEBUG("Finalize successful");
    return StatusCode::SUCCESS;
  }

  Amg::Vector3D TrackDensitySeedFinder::findSeed(const std::vector<const Trk::Track*> & VectorTrk,
						 const xAOD::Vertex * constraint) const
  {
    
    //create perigees from track list
    std::vector<const TrackParameters*> perigeeList;
    
    std::vector<const Trk::Track*>::const_iterator begin=VectorTrk.begin();
    std::vector<const Trk::Track*>::const_iterator end=VectorTrk.end();
    
    for (std::vector<const Trk::Track*>::const_iterator iter=begin;
	 iter!=end;++iter) 
    {
      perigeeList.push_back((*iter)->perigeeParameters());
    }
   
    //create seed from perigee list
    return findSeed(perigeeList,constraint);  
  }

  Amg::Vector3D TrackDensitySeedFinder::findSeed(const std::vector<const Trk::TrackParameters*> & perigeeList,
						 const xAOD::Vertex * constraint) const
  {
    double zResult {0.};
    if ( perigeeList.size()>0 ) 
    {
      zResult = m_densityEstimator->globalMaximum (perigeeList);
      ATH_MSG_DEBUG("Highest density Z position found: " << zResult);
    }
    else
    {
      ATH_MSG_DEBUG("No tracks with sufficient weight; return z position = 0");
    }

    if (constraint)
    {
      return Amg::Vector3D(constraint->position().x(), constraint->position().y(), zResult + constraint->position().z());
    }
    else
    {
      return Amg::Vector3D(0.,0.,zResult);
    }  
  }

    // testing find seed with width
  std::pair<Amg::Vector3D,Amg::MatrixX> TrackDensitySeedFinder::findAnalyticSeed(const std::vector<const Trk::TrackParameters*>& perigeeList, const xAOD::Vertex * constraint ) const{
      std::pair<double,double> zResult {0.,0.};
    if ( perigeeList.size()>0 ) 
    {
      zResult = m_densityEstimator->globalMaximumWithWidth (perigeeList);
      ATH_MSG_DEBUG("Highest density Z position found: " << zResult.first << "with width: " << zResult.second);
    }
    else
    {
      ATH_MSG_DEBUG("No tracks with sufficient weight; return z position = 0");
    }

    if (constraint)
    {
      Amg::Vector3D positionVector = Amg::Vector3D(constraint->position().x(), constraint->position().y(), zResult.first + constraint->position().z());
      Amg::MatrixX covarianceMatrix = constraint->covariancePosition() ;
      
      ATH_MSG_DEBUG("The vertex seed width is " << zResult.second);
     
     // TODO: not sure what to do when no seed is found as width is NaN...
      // if seed position is 0 then vertexing should stop anayways...
      if(zResult.second!=zResult.second) covarianceMatrix.fillSymmetric(2,2,1); 
      else covarianceMatrix.fillSymmetric(2,2,std::pow(zResult.second,2.));

      return std::make_pair(positionVector,covarianceMatrix) ;
    }
    else
    {
      Amg::Vector3D positionVector = Amg::Vector3D(0.,0.,zResult.first);
      Amg::MatrixX covarianceMatrix;
      for(const auto i : {0,1,2} ){
        for(const auto j : {0,1,2}){
          covarianceMatrix(i,j)=0;
        }
      }
      return {positionVector,covarianceMatrix};
    }  
}

  std::vector<Amg::Vector3D> TrackDensitySeedFinder::findMultiSeeds(const std::vector<const Trk::Track*>& /* vectorTrk */,
								    const xAOD::Vertex * /* constraint */)  const
  {
    //implemented to satisfy inheritance but this algorithm only supports one seed at a time
    ATH_MSG_WARNING("Multi-seeding requested but seed finder not able to operate in that mode, returning no seeds");
    return std::vector<Amg::Vector3D>(0);
  }

  std::vector<Amg::Vector3D> TrackDensitySeedFinder::findMultiSeeds(const std::vector<const Trk::TrackParameters*>& /* perigeeList */,
								    const xAOD::Vertex * /* constraint */) const
  {
     //implemented to satisfy inheritance but this algorithm only supports one seed at a time
    ATH_MSG_WARNING("Multi-seeding requested but seed finder not able to operate in that mode, returning no seeds");
    return std::vector<Amg::Vector3D>(0);
  }


}
