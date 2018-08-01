#pragma once

#include "Acts/Seeding/SeedSPGrid.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <vector>
#include <memory>
#include <queue>

namespace Acts{
     
    struct LinCircle{
      float Zo;
      float cotTheta;
      float iDeltaR;
      float Er;
      float U;
      float V;
    };

    struct SeedmakerState {
      // grid with ownership of all InternalSpacePoint
      std::unique_ptr<Acts::SPGrid> binnedSP;

      // contains parameters required to calculate circle with linear equation 
      // ...for bottom-middle
      std::vector<LinCircle> linCircleBottom;
      // ...for middle-top
      std::vector<LinCircle> linCircleTop;

      // middle spacepoint bin indices
      size_t phiIndex = 1;
      size_t zIndex = 1;

      // container with seeds created so far
      std::queue<std::unique_ptr<const InternalSeed> > outputQueue;
    };
}
