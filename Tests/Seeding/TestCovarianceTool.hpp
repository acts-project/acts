#pragma once


/// @class ICovarianceTool
/// Interface Class that returns covariance on SPForSeed locations in r and z

#include "Acts/Seeding/ICovarianceTool.hpp"
#include "SpacePoint.hpp"

#include <boost/type_erasure/any_cast.hpp>


namespace Acts {
class CovarianceTool : public ICovarianceTool
{
public:
  /// Virtual destructor
    virtual ~CovarianceTool() = default; 
    
    /// ICovarianceTool interface method
    /// returns squared errors in z and r for the passed SpacePoint
    /// @param sp is the SpacePoint fro which the covariance values will be retrieved
    /// @param zAlign is the alignment uncertainty in z. 
    /// it is going to be squared and added to covz.
    /// @param rAlign is the alignment uncertainty in r.
    /// it is going to be squared and added to covr.
    /// @param sigma is multiplied with the combined alignment and covariance errors
    std::array<float,2>
    getCovariances(const Acts::concept::AnySpacePoint<>* sp, 
                   float zAlign = 0, 
                   float rAlign = 0, 
                   float sigma=1);

};
    inline
    std::array<float,2>
    CovarianceTool::getCovariances(const Acts::concept::AnySpacePoint<>* sp, 
                                   float zAlign,
                                   float rAlign,
                                   float sigma)
    {
      std::array<float,2> cov;
      cov[0] = ((boost::type_erasure::any_cast<SpacePoint>(*sp)).covr + rAlign*rAlign) * sigma;
      cov[1] = ((boost::type_erasure::any_cast<SpacePoint>(*sp)).covz + zAlign*zAlign) * sigma;
      return cov;
    }
}
