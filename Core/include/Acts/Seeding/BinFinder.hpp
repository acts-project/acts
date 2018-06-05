#pragma once

#include "Acts/Seeding/IBinFinder.hpp"

#include <set>

namespace Acts {
class BinFinder : public IBinFinder
{
public:
/// Virtual destructor
  ~BinFinder() = default;

   virtual
   std::set<size_t>
   findBins(size_t phiBin, size_t zBin, std::unique_ptr<SPGrid>& binnedSP);

};
}
