
#pragma once

#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/ISeedFilter.hpp"
#include "ACTS/Seeding/IQualityTool.hpp"

namespace Acts{
namespace Seeding{
  class ATLSeedFilter : public SeedFilter{
    public: 
    ATLSeedFilter(SeedFilterConfig cfg,
               std::shared_ptr<IQualityTool> qualityTool);

    SeedFilter() = delete;
    ~SeedFilter();

    virtual
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
    filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                         std::shared_ptr<SPForSeed> middleSP,
                         std::vector<std::shared_ptr<SPForSeed>> topSpVec,
                         std::vector<float> invHelixRadiusVec,
                         std::vector<float> impactParametersVec,
                         float zOrigin) override;

    virtual
    std::vector<std::shared_ptr<InternalSeed> >
    filterSeeds_1SpFixed(std::vector<std::pair<float, std::shared_ptr<InternalSeed> > > seedsPerSpM) override;

    virtual
    std::vector<std::shared_ptr<Seed> >
    filterSeeds_byRegion(std::vector<std::shared_ptr<InternalSeed> > seedsPerRegion) override;

    private:
    const SeedFilterConfig m_cfg;
    const std::shared_ptr<IQualityTool> m_qualityTool;
  };
}
}
