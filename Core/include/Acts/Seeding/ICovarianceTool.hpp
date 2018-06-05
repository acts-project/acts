#pragma once

/// @class ICovarianceTool
/// Interface Class that sets covariance on SPForSeed locations in r and z

#include "Acts/Seeding/SPForSeed.hpp"

namespace Acts {
class ICovarianceTool
{
public:
  /// Virtual destructor
    virtual ~ICovarianceTool() = default; 
    
    /// ICovarianceTool interface method
    /// returns squared errors in z and r for the passed SpacePoint
    /// @param sp is the SpacePoint fro which the covariance values will be retrieved
    /// @param zAlign is the alignment uncertainty in z. 
    /// it is going to be squared and added to covz.
    /// @param rAlign is the alignment uncertainty in r.
    /// it is going to be squared and added to covr.
    /// @param sigma is multiplied with the combined alignment and covariance errors
    virtual std::array<float,2> getCovariances(const Acts::concept::AnySpacePoint<>* sp, float zAlign = 0, float rAlign = 0, float sigma=1) =0;

};
}
