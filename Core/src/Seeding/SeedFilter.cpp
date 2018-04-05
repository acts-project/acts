#include "ACTS/Seeding/SeedFilter.hpp"
#include <utility>

namespace Acts{
namespace Seeding{
  //constructor
  SeedFilter::SeedFilter(float minQuality) : m_minQuality(minQuality){}
  
  //destructor
  SeedFilter::~SeedFilter(){}

  // function to filter seeds based on all seeds with same bottom- and middle-spacepoint.
  // return vector must contain quality of each seed 
  std::vector<std::pair<float, InternalSeed> >
  SeedFilter::filterSeeds_2SpFixed(SPForSeed* bottomSP,
                                         SPForSeed* middleSP,
                                         std::vector<SPForSeed*> topSpVec,
                                         std::vector<float> invHelixRadiusVec,
                                         std::vector<float> impactParametersVec,
                                         float zOrigin,
                                         std::unique_ptr<Acts::Seeding::Config>& config){
  
    std::vector<std::pair<float, InternalSeed> > selectedSeeds;
  
    float startquality=0.;
    bool SCTSeed = bottomSP->spacepoint->clusterList().second;
    // SCT seeds have highest purity and get highest quality 
    if (SCTSeed) startquality  = 400.;
    // if two compatible seeds with high distance in r are found, compatible seeds span 5 layers
    // -> very good seed
    float foundCompatibleRadius;
    bool foundCompatibleSeed;
  
    
    for(int i = 0; i < topSpVec.size(); i++){

      foundCompatibleSeed = false;
      float invHelixRadius = invHelixRadiusVec.at(i);
      float lowerLimitCurv = invHelixRadius - config->deltaInvHelixRadius; 
      float upperLimitCurv = invHelixRadius + config->deltaInvHelixRadius; 
      float currentTop_r   = topSpVec.at(i)->radius();
      float impact         = impactParametersVec.at(i);
      // TODO: check if filtering worse if quality without factor for deltaTheta
      float quality = startquality - impact;
      // pixel seeds are still better than mixed seeds, add to quality
      bool pixelSeed = !topSpVec.at(i)->spacepoint->clusterList().second;
      if (pixelSeed) quality += 200.;
      for (int j = 0; j < topSpVec.size(); j++){
        if (i == j) continue;
        // compared top SP should have at least deltaRMin distance
        float otherTop_r = topSpVec.at(j)->radius();
        float deltaR = currentTop_r - otherTop_r;
        if (std::abs(deltaR) < config->deltaRMin)   continue;
        // curvature difference within limits?
        // TODO: how much slower than sorting all vectors by curvature
        // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
        if (invHelixRadiusVec.at(j) < lowerLimitCurv) continue;
        if (invHelixRadiusVec.at(j) > upperLimitCurv) continue;
        if(foundCompatibleSeed){
          if (std::abs(foundCompatibleRadius - otherTop_r) > config->deltaRMin){ quality +=200.; break;}
        } else {
          foundCompatibleSeed = true;
          foundCompatibleRadius = otherTop_r;
          quality+=200.;
        }
      }
      // discard low quality seeds
      if (quality < m_minQuality) continue;
      // mixed seeds have low purity, so only accept mixed seeds with at least one compatible seed
      // FIXME: test for space point quality removed. test how much more mixed seeds we have.
      if(pixelSeed == SCTSeed && quality < 0.) continue;
      selectedSeeds.push_back(std::make_pair(quality, InternalSeed(bottomSP,middleSP,topSpVec.at(i),zOrigin)));
      }
    return selectedSeeds;
  }

  // after creating all seeds with a common middle space point, filter again
  void SeedFilter::filterSeeds_1SpFixed(std::unique_ptr<Acts::Seeding::Cache>& cache,
                                        std::unique_ptr<Acts::Seeding::Config>& config){
  
    std::sort(cache->seedsPerSpM.begin(),cache->seedsPerSpM.end(),comQuality());
    int itEnd = cache->seedsPerSpM.size();
    if(itEnd > config->maxSeedsPerSpM){
      itEnd = config->maxSeedsPerSpM + 1;
    }
  
    typename std::vector<std::pair<float,std::shared_ptr<InternalSeed>> >::iterator
    lf = cache->seedsPerSpM.begin(),
    l  = cache->seedsPerSpM.begin(),
    le = cache->seedsPerSpM.begin()+itEnd;
  
    if(l==le) return;
  
    std::shared_ptr<InternalSeed> s;
  
    for(; l!=le; ++l) {
  //  Testing various quality measurements to find out if seed should be included in output collection
      float w = (*l).first ;
      s       = (*l).second;
  //  TODO: This is very ATLAS specific
  //  if not first seed, if bottom SP is on IBL and quality w worse than 200  then ignore seed
      if(l!=lf && s->spacepoint0()->radius() < 43. && w < 200.) continue;
  //  if seed is mixed and quality w is worse than all SP qualities, ignore seed.
  //  sets quality to spacepoints & to seed, return bool used to skip mixed, low-quality seeds
      bool SSS = s->spacepoint0()->spacepoint->clusterList().second;
      bool PPP = !s->spacepoint2()->spacepoint->clusterList().second;
      // skip mixed seed if all SPs have already been used in higher quality pure seed
      if(!SSS && !PPP){
        if (s->spacepoint0()->quality() > w && s->spacepoint1()->quality() > w && s->spacepoint2()->quality() > w ){
          continue;
        }
      }
      s->setQuality(w, PPP||SSS); 
  //  add dominating weight to all quality factors to order them by Seed type.
  //  used to sort output in this order: SSS, SSP, SPP, PPP
  //  in filterAllSeeds this leads to: mixed seeds are filtered if its quality
  //  not better than worst SSS seeds using these SP, but not PPP
      if     (s->spacepoint0()->spacepoint->clusterList().second) w+=3000.;
      else if(s->spacepoint1()->spacepoint->clusterList().second) w+=2000.;
      else if(s->spacepoint2()->spacepoint->clusterList().second) w+=1000.;
      cache->regionSeeds.push_back(std::make_pair(w,s));
    }
    cache->seedsPerSpM.clear();
  }

  void SeedFilter::filterSeeds_byRegion(std::unique_ptr<Acts::Seeding::Cache>& cache){
    
    std::sort(cache->regionSeeds.begin(), cache->regionSeeds.end(), comQuality()); 
    for(auto weightSeed : cache->regionSeeds){
      std::shared_ptr<InternalSeed> s = weightSeed.second;
      bool SSS = s->spacepoint0()->spacepoint->clusterList().second;
      bool PPP = !(s->spacepoint2()->spacepoint->clusterList().second);
      float w = weightSeed.first;
      if(!SSS && !PPP){
        if (s->spacepoint0()->quality() > w && s->spacepoint1()->quality() > w && s->spacepoint2()->quality() > w ){
          continue;
        }
      }
      s->setQuality(w, PPP||SSS);  
      Seed outSeed = Seed(s->spacepoint0()->spacepoint,
          s->spacepoint1()->spacepoint,
          s->spacepoint2()->spacepoint,
          s->quality());
      cache->outputSeeds.push_back(std::make_unique<Seed>(outSeed ) );
    }
    cache->regionSeeds.clear();
  }
}//namespace Seeding
}//namespace Acts
