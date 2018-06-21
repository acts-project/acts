// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    /// @param sp is the SpacePoint for which the covariance values will be
    /// retrieved
    /// @param zAlign is the alignment uncertainty in z. 
    /// it is going to be squared and added to covz.
    /// @param rAlign is the alignment uncertainty in r.
    /// it is going to be squared and added to covr.
    /// @param sigma is multiplied with the combined alignment and covariance
    /// errors
    /// @returns squared errors in z and r for the passed SpacePoint (ignoring
    /// correlations)
    virtual Acts::Vector2D getCovariances(const Acts::concept::AnySpacePoint<>* sp, float zAlign = 0, float rAlign = 0, float sigma=1) =0;

};
}
