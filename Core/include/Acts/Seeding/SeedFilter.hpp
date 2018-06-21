// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <memory>
#include <queue>

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/IQualityTool.hpp"
#include "Acts/Seeding/Seed.hpp"

namespace Acts{
  struct SeedFilterConfig{
    // the allowed delta between two inverted seed radii for them to be considered compatible.
    float deltaInvHelixRadius = 0.00003;
    // the impact parameters (d0) is multiplied by this factor and subtracted from quality
    float impactQualityFactor = 1.;
    // seed quality increased by this value if a compatible seed has been found.
    float compatSeedQuality = 200.;
    // minimum distance between compatible seeds to be considered for quality boost
    float deltaRMin = 5.;
    // in dense environments many seeds may be found per middle space point.
    // only seeds with the highest quality will be kept if this limit is reached.
    int maxSeedsPerSpM = 10;
    // how often do you want to increase the quality of a seed for finding a compatible seed?
    int compatSeedLimit = 2;
  };

/// @class Filter seeds at various stages with the currently
/// available information.
  class SeedFilter{
    public: 
    SeedFilter(SeedFilterConfig cfg,
               std::shared_ptr<IQualityTool> qualityTool);

    SeedFilter() = delete;
    ~SeedFilter();

/// Create InternalSeeds for the all seeds with the same bottom and middle
/// space point and discard all others.
/// @param bottomSP fixed bottom space point
/// @param middleSP fixed middle space point
/// @param topSpVec vector containing all space points that may be compatible
/// with both bottom and middle space point
/// @param origin on the z axis as defined by bottom and middle space point
/// @return vector of pairs containing seed quality and seed for all valid
/// created seeds
    virtual
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
    filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                         std::shared_ptr<SPForSeed> middleSP,
                         std::vector<std::shared_ptr<SPForSeed>>& topSpVec,
                         std::vector<float>& invHelixRadiusVec,
                         std::vector<float>& impactParametersVec,
                         float zOrigin) const;

/// Filter seeds once all seeds for one middle space point have been created
/// @param seedsPerSpM vector of pairs containing quality and seed for all
/// for all seeds with the same middle space point
/// @return vector of all InternalSeeds that not filtered out
    virtual
    void
    filterSeeds_1SpFixed(std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >& seedsPerSpM, std::queue<std::shared_ptr<InternalSeed> >& queue) const;

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
