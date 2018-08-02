// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "Acts/Seeding/SeedFilter.hpp"
#include <utility>

namespace Acts{
  //constructor
  SeedFilter::SeedFilter(SeedFilterConfig config,
                         IExperimentCuts* expCuts)
                         :m_cfg (config),
                          m_experimentCuts (expCuts){}

  SeedFilter::SeedFilter(SeedFilterConfig config): m_cfg (config) {}

  //destructor
  SeedFilter::~SeedFilter(){}

  // function to filter seeds based on all seeds with same bottom- and middle-spacepoint.
  // return vector must contain weight of each seed 
  std::vector<std::pair<float, std::unique_ptr<const InternalSeed> > >
  SeedFilter::filterSeeds_2SpFixed(const InternalSpacePoint* bottomSP,
                                   const InternalSpacePoint* middleSP,
                                   std::vector<const InternalSpacePoint* >& topSpVec,
                                   std::vector<float>& invHelixRadiusVec,
                                   std::vector<float>& impactParametersVec,
                                   float zOrigin) const {
  
    std::vector<std::pair<float, std::unique_ptr<const InternalSeed> > > selectedSeeds;
  
    // if two compatible seeds with high distance in r are found, compatible seeds span 5 layers
    // -> very good seed
    std::vector<float> compatibleSeedRadii;

    for(size_t i = 0; i < topSpVec.size(); i++){

      float invHelixRadius = invHelixRadiusVec[i];
      float lowerLimitCurv = invHelixRadius - m_cfg.deltaInvHelixRadius; 
      float upperLimitCurv = invHelixRadius + m_cfg.deltaInvHelixRadius; 
      float currentTop_r   = topSpVec[i]->radius();
      float impact         = impactParametersVec[i];

      float weight = -(impact * m_cfg.impactWeightFactor);
      for (size_t j = 0; j < topSpVec.size(); j++){
        if (i == j) continue;
        // compared top SP should have at least deltaRMin distance
        float otherTop_r = topSpVec[j]->radius();
        float deltaR = currentTop_r - otherTop_r;
        if (std::abs(deltaR) < m_cfg.deltaRMin)   continue;
        // curvature difference within limits?
        // TODO: how much slower than sorting all vectors by curvature
        // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
        if (invHelixRadiusVec[j] < lowerLimitCurv) continue;
        if (invHelixRadiusVec[j] > upperLimitCurv) continue;
        bool newCompSeed = true;
        for(float previousRadius : compatibleSeedRadii){
          // original ATLAS code uses higher min distance for 2nd found compatible seed (20mm instead of 5mm)
          if(std::abs(previousRadius - otherTop_r) < m_cfg.deltaRMin) {newCompSeed = false; break;}
        }
        if(newCompSeed)
        {
          compatibleSeedRadii.push_back(otherTop_r);
          weight+= m_cfg.compatSeedWeight;
        }
        if(compatibleSeedRadii.size() >= m_cfg.compatSeedLimit) break;
      }
      if(m_experimentCuts != nullptr){
        // add detector specific considerations on the seed weight
        weight += m_experimentCuts->seedWeight(bottomSP, middleSP, topSpVec[i]);
        // discard seeds according to detector specific cuts (e.g.: weight)
        if (!m_experimentCuts->singleSeedCut(weight, bottomSP, middleSP, topSpVec[i])) continue;
      }
      selectedSeeds.push_back(std::make_pair(weight, std::make_unique<const InternalSeed>(bottomSP,middleSP,topSpVec[i],zOrigin)));
      }
    return selectedSeeds;
  }



  // after creating all seeds with a common middle space point, filter again
  void
  SeedFilter::filterSeeds_1SpFixed(std::vector<std::pair<float,std::unique_ptr<const InternalSeed > > >& seedsPerSpM, std::queue<std::unique_ptr<const InternalSeed> >& queue, std::mutex& outputMutex) const {

    //sort by weight and iterate only up to configured max number of seeds per middle SP
    std::sort((seedsPerSpM.begin()),(seedsPerSpM.end()),[]
                   (const std::pair<float,std::unique_ptr<const Acts::InternalSeed>>& i1,
                    const std::pair<float,std::unique_ptr<const Acts::InternalSeed>>& i2)
                   {
                     return i1.first > i2.first;
                   });
    if(m_experimentCuts != nullptr){
      seedsPerSpM = m_experimentCuts->cutPerMiddleSP(std::move(seedsPerSpM));
    }
    unsigned int maxSeeds = seedsPerSpM.size();
    if(maxSeeds > m_cfg.maxSeedsPerSpM){
      maxSeeds = m_cfg.maxSeedsPerSpM + 1;
    }
    auto itBegin = seedsPerSpM.begin();
    auto it = seedsPerSpM.begin();
    // default filter removes the last seeds if maximum amount exceeded
    // ordering by weight by filterSeeds_2SpFixed means these are the lowest weight seeds
    outputMutex.lock();
    for(; it<itBegin+maxSeeds; ++it) {
      queue.push(std::move((*it).second));
    }
    outputMutex.unlock();
  }

}//namespace Acts
