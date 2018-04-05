#ifndef ACTS_SEEDING_ATL_COVARIANCETOOL_H
#define ACTS_SEEDING_ATL_COVARIANCETOOL_H
#include "ACTS/Seeding/ICovarianceTool.hpp"

namespace Acts {

class SPForSeed;

/// @class ATL_CovarianceTool
///
/// ICovarianceTool implementation that sets covariance
/// on SPForSeed locations in r and z and expects
/// SpacePoint to be ATLAS
/// SpacePoints to have a covariance() method, which returns a 
/// covariance matrix with error in r (endcaps) resp. z (barrel) at index 1,1
class ATL_CovarianceTool : public ICovarianceTool 
{
public:

  /// Constructor
  virtual ATL_CovarianceTool() = default;
  /// Destructor
  virtual ~ATL_CovarianceTool() = default;

  /// sets squared errors in z and r for the passed SPForSeed
  /// @param sp is the SPForSeed whose covariance values will be set
  /// @param zAlign is the alignment uncertainty in z. 
  /// it is going to be squared and added to covz.
  /// @param rAlign is the alignment uncertainty in r.
  /// it is going to be squared and added to covr.
  void setCovariances(SPForSeed& sp, float zAlign = 0, float rAlign = 0,float sigmaError=1) final;
}

#include "ACTS/Seeding/ATL_CovarianceTool.ipp"
#endif  // ACTS_SEEDING_COVARIANCETOOL_H
