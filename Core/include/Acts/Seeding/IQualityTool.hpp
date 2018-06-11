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
/// @class IQualityTool can be used to increase or decrease seed weights
/// (quality) based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed quality cut. As the quality is stored in
/// seeds, there are two distinct methods.
    class IQualityTool{
    public:
/// Returns seed quality bonus/malus depending on detector considerations.
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return seed weight to be added to the seed's quality
      virtual float seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) = 0;
/// @param quality the current seed quality
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return true if the seed should be kept, false if the seed should be
/// discarded
      virtual bool passesQualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) = 0;
    };
}
