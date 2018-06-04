#ifndef IBINFINDER_H
#define IBINFINDER_H

#include "ACTS/Seeding/SeedSPGrid.hpp"


namespace Acts {
class IBinFinder
{
public:
  /// Virtual destructor
    virtual
    ~IBinFinder() = default;

    virtual
    std::set<size_t>
    findBins(size_t phiBin,size_t zBin, std::unique_ptr<SPGrid>& binnedSP) =0;

};
}
#endif //IBINFINDER_H

