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
      std::unique_ptr<Acts::SPGrid> binnedSP;
      std::vector<InternalSeed> seeds;

      std::vector<LinCircle> linCircleBottom;
      std::vector<LinCircle> linCircleTop;

      float highland;
      float maxScatteringAngle2;
      float pTPerHelixRadius;
      float minHelixRadius2;


      size_t phiIndex = 1;
      size_t zIndex = 1;

      std::queue<std::shared_ptr<Seed>> outputQueue;


    };
}
