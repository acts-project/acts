#ifndef ISEEDINGBIN_H
#define ISEEDINGBIN_H

#include "ACTS/Seeding/ISeedingBin.hpp"

namespace Acts {
namespace Seeding {
class ISeedingBin
{
public:
  /// Virtual destructor
    virtual ~ISeedingBin() = default;

    virtual void findBins(ISeedingBin bin) =0;

};
}
}
#endif //ISEEDINGBIN_H
