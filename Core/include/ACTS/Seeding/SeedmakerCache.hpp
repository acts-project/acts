#ifndef SeedmakerCache_hpp
#define SeedmakerCache_hpp

#include "ACTS/Seeding/SeedSPGrid.hpp"
#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/Seed.hpp"

#include <vector>
#include <memory>

namespace Acts{
     
    struct LinCircle{
      float Zo;
      float cotTheta;
      float iDeltaR;
      float Er;
      float U;
      float V;
    };

    struct Cache {
      std::unique_ptr<Acts::SPGrid> binnedSP;
      std::vector<InternalSeed> seeds;


      std::vector<LinCircle> linCircleBottom;
      std::vector<LinCircle> linCircleTop;
      std::vector<std::shared_ptr<Seed> > outputSeeds;

      float highland;
      float maxScatteringAngle2;
      float pTPerHelixRadius;
      float minHelixRadius2;


    };
}
#endif //SeedmakerCache_hpp
