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
                         std::shared_ptr<IQualityTool> qualityTool)
                         :m_cfg (config),
                          m_qualityTool (qualityTool){}

  SeedFilter::SeedFilter(SeedFilterConfig config): m_cfg (config) {}

  //destructor
  SeedFilter::~SeedFilter(){}

  // function to filter seeds based on all seeds with same bottom- and middle-spacepoint.
  // return vector must contain quality of each seed 
  // TODO: performance of shared_ptr vs raw pointer?
  std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
  SeedFilter::filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                                   std::shared_ptr<SPForSeed> middleSP,
                                   std::vector<std::shared_ptr<SPForSeed >>& topSpVec,
                                   std::vector<float>& invHelixRadiusVec,
                                   std::vector<float>& impactParametersVec,
                                   float zOrigin) const {
  
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > > selectedSeeds;
  
    // if two compatible seeds with high distance in r are found, compatible seeds span 5 layers
    // -> very good seed
    std::vector<float> compatibleSeedRadii;

    for(size_t i = 0; i < topSpVec.size(); i++){

      float invHelixRadius = invHelixRadiusVec.at(i);
      float lowerLimitCurv = invHelixRadius - m_cfg.deltaInvHelixRadius; 
      float upperLimitCurv = invHelixRadius + m_cfg.deltaInvHelixRadius; 
      float currentTop_r   = topSpVec.at(i)->radius();
      float impact         = impactParametersVec.at(i);

      float quality = -(impact * m_cfg.impactQualityFactor);
      for (size_t j = 0; j < topSpVec.size(); j++){
        if (i == j) continue;
        // compared top SP should have at least deltaRMin distance
        float otherTop_r = topSpVec.at(j)->radius();
        float deltaR = currentTop_r - otherTop_r;
        if (std::abs(deltaR) < m_cfg.deltaRMin)   continue;
        // curvature difference within limits?
        // TODO: how much slower than sorting all vectors by curvature
        // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
        if (invHelixRadiusVec.at(j) < lowerLimitCurv) continue;
        if (invHelixRadiusVec.at(j) > upperLimitCurv) continue;
        bool newCompSeed = true;
        for(float previousRadius : compatibleSeedRadii){
          // original ATLAS code uses higher min distance for 2nd found compatible seed (20mm instead of 5mm)
          if(std::abs(previousRadius - otherTop_r) < m_cfg.deltaRMin) {newCompSeed = false; break;}
        }
        if(newCompSeed)
        {
          compatibleSeedRadii.push_back(otherTop_r);
          quality+= m_cfg.compatSeedQuality;
        }
        if(compatibleSeedRadii.size() >= m_cfg.compatSeedLimit) break;
      }
      if(m_qualityTool != nullptr){
        // add detector specific considerations on the seed quality
        quality += m_qualityTool->seedQuality(bottomSP, middleSP, topSpVec.at(i));
        // discard low quality seeds
        if (!m_qualityTool->passesQualityCut(quality, bottomSP, middleSP, topSpVec.at(i))) continue;
      }
      selectedSeeds.push_back(std::make_pair(quality, std::make_shared<InternalSeed>(bottomSP,middleSP,topSpVec.at(i),zOrigin)));
      }
    return selectedSeeds;
  }



  // after creating all seeds with a common middle space point, filter again
  // TODO: performance of shared_ptr vs raw pointer?
  void
  SeedFilter::filterSeeds_1SpFixed(std::vector<std::pair<float,std::shared_ptr<InternalSeed > > >& seedsPerSpM, std::queue<std::shared_ptr<InternalSeed> >& queue) const {

    //sort by quality and iterate only up to configured max number of seeds per middle SP
    std::sort(seedsPerSpM.begin(),seedsPerSpM.end(),comQuality());
    unsigned int maxSeeds = seedsPerSpM.size();
    if(maxSeeds > m_cfg.maxSeedsPerSpM){
      maxSeeds = m_cfg.maxSeedsPerSpM + 1;
    }
    auto itBegin = seedsPerSpM.begin();
    auto it = seedsPerSpM.begin();
  
    std::vector<std::shared_ptr<InternalSeed> > filteredSeeds;
    // default filter removes the last seeds if maximum amount exceeded
    // ordering by quality by filterSeeds_2SpFixed means these are the lowest quality seeds
    for(; it<itBegin+maxSeeds; ++it) {
      std::shared_ptr<InternalSeed> internalSeed       = (*it).second;

      //auto outSeed = std::make_shared<Seed>(internalSeed->spacepoint0()->spacepoint,
      //                                      internalSeed->spacepoint1()->spacepoint,
      //                                      internalSeed->spacepoint2()->spacepoint,
      //                                      internalSeed->z());
      queue.push(internalSeed);
    }
  }

}//namespace Acts
