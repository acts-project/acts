#include "ACTS/Seeding/BinFinder.hpp"

//std::set<std::array<int, 2> >
std::set<size_t>
Acts::Seeding::BinFinder::findBins(size_t i, size_t j, std::unique_ptr<Acts::Seeding::SPGrid>& binnedSP){
  return binnedSP->neighborHoodIndices({i,j});  
}
