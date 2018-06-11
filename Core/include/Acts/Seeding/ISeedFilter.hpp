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

namespace Acts {
/// @class Interface to filter seeds at various stages with the currently
/// available information.
  class ISeedFilter{
    public: 
/// default destructor
    virtual
    ~ISeedFilter(){}
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
                         std::vector<std::shared_ptr<SPForSeed> >& topSpVec,
                         std::vector<float>& invHelixRadiusVec,
                         std::vector<float>& impactParametersVec,
                         float zOrigin) const = 0;
/// Filter seeds once all seeds for one middle space point have been created
/// @param seedsPerSpM vector of pairs containing quality and seed for all
/// for all seeds with the same middle space point
/// @return vector of all InternalSeeds that not filtered out
    virtual
    void
    filterSeeds_1SpFixed(std::vector<std::pair<float,std::shared_ptr<InternalSeed > > >& seedsPerSpM, std::queue<std::shared_ptr<Seed> >& queue) const = 0;
  };
}
