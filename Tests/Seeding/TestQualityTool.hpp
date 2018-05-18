#pragma once

#include "ACTS/Seeding/IQualityTool.hpp"
#include <memory>

namespace Acts{
  namespace Seeding{
    class QualityTool : public IQualityTool{
    public:
      float seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) override;
      bool passesQualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top) override;
    };
    inline float QualityTool::seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top){
      return 0.;
    } 
    inline bool QualityTool::passesQualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top){
      return true;
    }
  }
}
