#ifndef SeedmakerCache_hpp
#define SeedmakerCache_hpp

#include "ACTS/Seeding/SeedSPGrid.hpp"
#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/Seed.hpp"

#include <vector>
#include <memory>

namespace Acts{
  namespace Seeding{
     
    struct LinCircle{
      float Zo;
      float cotTheta;
      float iDeltaR;
      float Er;
      float U;
      float V;
    };

    struct Cache {
      std::unique_ptr<Acts::Seeding::SPGrid> binnedSP;
      std::vector<InternalSeed> seeds;

      std::vector<std::pair<float,std::shared_ptr<InternalSeed > > > seedsPerSpM;
      std::vector<std::pair<float,std::shared_ptr<InternalSeed > > > regionSeeds;

      std::vector<LinCircle> linCircleBottom;
      std::vector<LinCircle> linCircleTop;
      std::vector<std::unique_ptr<Seed> > outputSeeds;

      float highland;
      float maxScatteringAngle2;
      float pTPerHelixRadius;
      float minHelixRadius2;


    };
  }
}
#endif //SeedmakerCache_hpp
