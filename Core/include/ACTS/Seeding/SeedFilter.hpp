#pragma once

#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/ISeedFilter.hpp"
#include "ACTS/Seeding/IQualityTool.hpp"

namespace Acts{
  struct SeedFilterConfig{
    // the allowed delta between two seed radii for them to be considered compatible.
    float deltaInvHelixRadius = 0.00003;
    // the impact parameters (d0) is multiplied by this factor and subtracted from quality
    float impactQualityFactor = 1.;
    // seed quality increased by this value if a compatible seed has been found.
    float compatSeedQuality = 200.;
    // minimum distance between compatible seeds to be considered for quality boost
    float deltaRMin = 5.;
    // in jets many seeds may be found per middle space point.
    // only seeds with the highest quality will be kept if this limit is reached.
    int maxSeedsPerSpM = 10;
    // how often do you want to increase the quality of a seed for finding a compatible seed?
    int compatSeedLimit = 2;
  };

  class SeedFilter : public ISeedFilter{
    public: 
    SeedFilter(SeedFilterConfig cfg,
               std::shared_ptr<IQualityTool> qualityTool);

    SeedFilter() = delete;
    ~SeedFilter();

    virtual
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
    filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                         std::shared_ptr<SPForSeed> middleSP,
                         std::vector<std::shared_ptr<SPForSeed>>& topSpVec,
                         std::vector<float>& invHelixRadiusVec,
                         std::vector<float>& impactParametersVec,
                         float zOrigin) const override;

    virtual
    std::vector<std::shared_ptr<InternalSeed> >
    filterSeeds_1SpFixed(std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >& seedsPerSpM) const override;

    virtual
    std::vector<std::shared_ptr<Seed> >
    filterSeeds_byRegion(std::vector<std::shared_ptr<InternalSeed> >& seedsPerRegion) const override;

    private:
    const SeedFilterConfig m_cfg;
    const std::shared_ptr<IQualityTool> m_qualityTool;
  };
}


// quality comparator, returns true if first value is larger than second value
// using comQuality to sort will return higher quality -> lower quality.
class comQuality  {
public:
  
  bool operator ()
  (const std::pair<float,std::shared_ptr<Acts::InternalSeed>>& i1,
   const std::pair<float,std::shared_ptr<Acts::InternalSeed>>& i2)
  {
    return i1.first > i2.first;
  }
};
