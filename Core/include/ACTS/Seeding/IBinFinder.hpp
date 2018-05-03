#ifndef IBINFINDER_H
#define IBINFINDER_H

#include "ACTS/Seeding/SeedSPGrid.hpp"


namespace Acts {
namespace Seeding{
class IBinFinder
{
public:
  /// Virtual destructor
    virtual
    ~IBinFinder() = default;

    virtual
    std::set<size_t>
    findBins(int phiBin,int zBin, std::unique_ptr<SPGrid>& binnedSP) =0;

};
}
}
#endif //IBINFINDER_H

