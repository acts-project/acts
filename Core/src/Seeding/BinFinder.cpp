#include "Acts/Seeding/BinFinder.hpp"

std::set<size_t>
Acts::BinFinder::findBins(size_t i, size_t j, std::unique_ptr<Acts::SPGrid>& binnedSP){
  return binnedSP->neighborHoodIndices({i,j});  
}
