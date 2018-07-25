#include "Acts/Seeding/IExperimentCuts.hpp"
#include <algorithm>


namespace Acts{
    class ATLASCuts : public IExperimentCuts{
    public:
/// Returns seed weight bonus/malus depending on detector considerations.
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return seed weight to be added to the seed's weight
      float seedWeight(const InternalSpacePoint* bottom, const InternalSpacePoint* middle, const InternalSpacePoint* top) const;
/// @param weight the current seed weight
/// @param bottom bottom space point of the current seed
/// @param middle middle space point of the current seed
/// @param top top space point of the current seed
/// @return true if the seed should be kept, false if the seed should be
/// discarded
      bool singleSeedCut(float weight, const InternalSpacePoint* bottom, const InternalSpacePoint* middle, const InternalSpacePoint* top) const;

/// @param seeds contains pairs of weight and seed created for one middle space
/// point
/// @return vector of seeds that pass the cut
      std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > > cutPerMiddleSP(std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > > seeds) const;
};

float
ATLASCuts::seedWeight(const InternalSpacePoint* bottom, const InternalSpacePoint*, const InternalSpacePoint* top) const {
  float weight = 0;
  if( bottom->radius() > 150){ weight = 400;}
  if( top->radius() < 150){ weight = 200;}
  return weight;
}
bool ATLASCuts::singleSeedCut(float, const InternalSpacePoint*, const InternalSpacePoint*, const InternalSpacePoint*) const {
  return true;
}

std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > >
ATLASCuts::cutPerMiddleSP(std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > > seeds) const {
  std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > > newSeedsVector;
  if (seeds.size()>1){
    newSeedsVector.push_back(std::move(seeds[0]));
    size_t itLength = std::min(seeds.size(),size_t(5));
    for(int i = 1; i < itLength; i++){
      if (seeds[i].first > 200. || seeds[i].second->spacepoint0()->radius() > 43.){
        newSeedsVector.push_back(std::move(seeds[i]));
      }
    }
    return std::move(newSeedsVector);
  }
  return std::move(seeds);
}
}
