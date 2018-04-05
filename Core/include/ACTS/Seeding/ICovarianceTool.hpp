#ifndef ICOVARIANCETOOL_H
#define ICOVARIANCETOOL_H

/// @class ICovarianceTool
/// Interface Class that sets covariance on SPForSeed locations in r and z

#include "ACTS/Seeding/SPForSeed.hpp"

namespace Acts {
namespace Seeding{
class ICovarianceTool
{
public:
  /// Virtual destructor
    virtual ~ICovarianceTool() = default; 
    
    /// ICovarianceTool interface method
    /// sets squared errors in z and r for the passed SPForSeed
    /// @param sp is the SPForSeed whose covariance values will be set
    /// @param zAlign is the alignment uncertainty in z. 
    /// it is going to be squared and added to covz.
    /// @param rAlign is the alignment uncertainty in r.
    /// it is going to be squared and added to covr.
    virtual void setCovariances(SPForSeed& sp, float zAlign = 0, float rAlign = 0, float sigma=1) =0;

};
}
}
#endif //ICOVARIANCETOOL_H
