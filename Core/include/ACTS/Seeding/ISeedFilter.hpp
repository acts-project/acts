#ifndef ISeedFilter_h
#define ISeedFilter_h

#include <vector>
#include <memory>

#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/SeedmakerCache.hpp"

namespace Acts {
namespace Seeding {
  class ISeedFilter{
    public: 
    virtual
    ~ISeedFilter(){}

  
    virtual
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
    filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                         std::shared_ptr<SPForSeed> middleSP,
                         std::vector<std::shared_ptr<SPForSeed> > topSpVec,
                         std::vector<float> invHelixRadiusVec,
                         std::vector<float> impactParametersVec,
                         float zOrigin) = 0;

    virtual
    std::vector<std::shared_ptr<InternalSeed> > 
    filterSeeds_1SpFixed(std::vector<std::pair<float,std::shared_ptr<InternalSeed > > > seedsPerSpM) = 0;

    virtual 
    std::vector<std::shared_ptr<Seed> >
    filterSeeds_byRegion(std::vector<std::shared_ptr<InternalSeed> > seedsPerRegion) = 0;
    
  };
}
}
#endif  // ISeedFilter_h
