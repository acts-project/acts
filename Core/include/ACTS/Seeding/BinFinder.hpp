#ifndef BINFINDER_H
#define BINFINDER_H

#include "ACTS/Seeding/IBinFinder.hpp"

#include <set>

namespace Acts {
namespace Seeding {
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
}
#endif //BOTTOMBINFINDER_H
