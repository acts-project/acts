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
    void 
    filterSeeds_1SpFixed(std::unique_ptr<Acts::Seeding::Cache>& cache,
                         std::unique_ptr<Acts::Seeding::Config>& config) = 0;
  
    virtual
    std::vector<std::pair<float, InternalSeed> >
    filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                         std::shared_ptr<SPForSeed> middleSP,
                         std::vector<std::shared_ptr<SPForSeed> > topSpVec,
                         std::vector<float> invHelixRadiusVec,
                         std::vector<float> impactParametersVec,
                         float zOrigin,
                         std::unique_ptr<Acts::Seeding::Config>& config) = 0;

    virtual 
    void 
    filterSeeds_byRegion(std::unique_ptr<Acts::Seeding::Cache>& cache) = 0;
    
  };
}
}
#endif  // ISeedFilter_h
