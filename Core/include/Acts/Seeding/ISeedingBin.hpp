#pragma once

#include "Acts/Seeding/ISeedingBin.hpp"

namespace Acts {
class ISeedingBin
{
public:
  /// Virtual destructor
    virtual ~ISeedingBin() = default;

    virtual void findBins(ISeedingBin bin) =0;

};
}
