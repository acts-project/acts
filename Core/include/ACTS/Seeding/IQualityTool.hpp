#pragma once

#include "ACTS/Seeding/InternalSeed.hpp"
#include <memory>

namespace Acts{
    class IQualityTool{
    public:
      virtual float seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) = 0;
      virtual bool passesQualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) = 0;
    };
}
