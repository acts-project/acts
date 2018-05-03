#include "ACTS/Seeding/ATLQualityTool.hpp"

#include "TrkSpacePoint/SpacePoint.h"

#include <boost/any.hpp>

namespace Acts{
  namespace Seeding{

    float 
    ATLQualityTool::seedQuality(std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top)
    {
      Trk::SpacePoint* bottomAtlSp = boost::any_cast(bottom->spacepoint);
      Trk::SpacePoint* topAtlSp = boost::any_cast(top->spacepoint);

      if (bottomAtlasSp->clusterList().second) return 400.;
      if (!topAtlasSp->clusterList().second) return 200.;
      return 0;
    }

    bool
    ATLQualityTool::qualityCut(float quality, std::shared_ptr<SPForSeed> bottom, std::shared_ptr<SPForSeed> middle, std::shared_ptr<SPForSeed> top){
      // tests for clusters only if first check passed. only works if minQuality is higher than mixedMinQuality 
      return quality > m_minQuality;
    }
  }
}
