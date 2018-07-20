// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include <memory>

namespace Acts{
/// @class IWeightTool can be used to increase or decrease seed weights
/// (weight) based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed weight cut. As the weight is stored in
/// seeds, there are two distinct methods.
    class IExperimentCuts{
    public:
/// Returns seed weight bonus/malus depending on detector considerations.
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return seed weight to be added to the seed's weight
      virtual float seedWeight(std::shared_ptr<InternalSpacePoint> bottom, std::shared_ptr<InternalSpacePoint> middle, std::shared_ptr<InternalSpacePoint> top) const = 0;
/// @param weight the current seed weight
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return true if the seed should be kept, false if the seed should be
/// discarded
      virtual bool singleSeedCut(float weight, std::shared_ptr<InternalSpacePoint> bottom, std::shared_ptr<InternalSpacePoint> middle, std::shared_ptr<InternalSpacePoint> top) const = 0;

/// @param seeds contains pairs of weight and seed created for one middle space
/// point
/// @return vector of seeds that pass the cut
      virtual std::vector<std::pair<float,std::shared_ptr<InternalSeed> > > cutPerMiddleSP(std::vector<std::pair<float,std::shared_ptr<InternalSeed> > > seeds) const = 0;
    };
}
