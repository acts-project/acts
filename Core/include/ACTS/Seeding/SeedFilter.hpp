#ifndef SeedFilter_h
#define SeedFilter_h

#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/ISeedFilter.hpp"

namespace Acts{
namespace Seeding{
  class SeedFilter : ISeedFilter{
    public: 
    SeedFilter(float minQuality);
    SeedFilter() = delete;
    ~SeedFilter();

    void 
    filterSeeds_1SpFixed(std::unique_ptr<Acts::Seeding::Cache>& cache,
                         std::unique_ptr<Acts::Seeding::Config>& config);

    std::vector<std::pair<float, InternalSeed> >
    filterSeeds_2SpFixed(SPForSeed* bottomSP,
                         SPForSeed* middleSP,
                         std::vector<SPForSeed*> topSpVec,
                         std::vector<float> invHelixRadiusVec,
                         std::vector<float> impactParametersVec,
                         float zOrigin,
                         std::unique_ptr<Acts::Seeding::Config>& config);


    void 
    filterSeeds_byRegion(std::unique_ptr<Acts::Seeding::Cache>& cache);
    private:

    float m_minQuality;
  };
}
}

// quality comparator, returns true if first value is larger than second value
// using comQuality to sort will return higher quality -> lower quality.
class comQuality  {
public:
  
  bool operator ()
  (const std::pair<float,std::shared_ptr<Acts::Seeding::InternalSeed>>& i1,
   const std::pair<float,std::shared_ptr<Acts::Seeding::InternalSeed>>& i2)
  {
    return i1.first > i2.first;
  }
};

#endif  // SeedFilter_h
