#pragma once

#include "ACTS/Seeding/IQualityTool.hpp"

namespace Acts{
  namespace Seeding{
    class ATLQualityTool : public IQualityTool{
    public:
      float seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) override;
      bool qualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) override;
    };
  }
}
