#ifndef BINFINDER_H
#define BINFINDER_H

#include "ACTS/Seeding/IBinFinder.hpp"

#include <set>

namespace Acts {
namespace Seeding {
class BinFinder : IBinFinder
{
public:
/// Virtual destructor
  ~BinFinder() = default;

//   std::set<std::array<int, 2> >
   std::set<size_t>
   findBins(int phiBin, int zBin, std::unique_ptr<SPGrid>& binnedSP);

};
}
}
#endif //BOTTOMBINFINDER_H
