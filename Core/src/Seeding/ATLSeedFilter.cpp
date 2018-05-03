#include "ACTS/Seeding/SeedFilter.hpp"
#include <utility>

namespace Acts{
namespace Seeding{
  ATLSeedFilter::ATLSeedFilter(SeedFilterConfig config,
                         std::shared_ptr<IQualityTool> qualityTool)
                         :SeedFilter(config, qualityTool){}



  // after creating all seeds with a common middle space point, filter again
  std::vector<std::pair<float,std::shared_ptr<InternalSeed > > >
  ATLSeedFilter::filterSeeds_1SpFixed(std::vector<std::pair<float,std::shared_ptr<InternalSeed > > > seedsPerSpM){

    //sort by quality and iterate only up to configured max number of seeds per middle SP
    int itEnd = seedsPerSpM.size();
    if(itEnd > m_cfg.maxSeedsPerSpM){
      std::sort(seedsPerSpM.begin(),seedsPerSpM.end(),comQuality());
      itEnd = m_cfg.maxSeedsPerSpM + 1;
    }
  
    typename std::vector<std::pair<float,std::shared_ptr<InternalSeed>> >::iterator
    lf = seedsPerSpM.begin(),
    l  = seedsPerSpM.begin(),
    le = seedsPerSpM.begin()+itEnd;
  
    std::vector<std::shared_ptr<InternalSeed> > filteredSeeds;
  
    for(; l!=le; ++l) {
  //  Testing various quality measurements to find out if seed should be included in output collection
      float w = (*l).first ;
      std::shared_ptr<InternalSeed> internalSeed       = (*l).second;
  //  TODO: this is the only line that is ATLAS specific
  //  if not first seed, if bottom SP is on IBL and quality w worse than 200  then ignore seed
      if(l!=lf && internalSeed->spacepoint0()->radius() < 43. && w < 200.) continue;
      filteredSeeds.push_back(internalSeed);
    }
    return filteredSeeds;
  }

} // namespace Seeding
} // namespace Acts
