// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {
class CovarianceTool {
 public:
  /// Virtual destructor
  ~CovarianceTool() = default;

  /// ICovarianceTool interface method
  /// returns squared errors in z and r for the passed SpacePoint
  /// @param sp is the SpacePoint from which the covariance values will be
  /// retrieved
  /// @param zAlign is the alignment uncertainty in z.
  /// it is going to be squared and added to varianceZ.
  /// @param rAlign is the alignment uncertainty in r.
  /// it is going to be squared and added to varianceR.
  /// @param sigma is multiplied with the combined alignment and covariance
  /// errors
  template <typename SpacePoint>
  Acts::Vector2 getCovariances(const SpacePoint* sp, float zAlign = 0,
                               float rAlign = 0, float sigma = 1);
};
template <typename SpacePoint>
inline Acts::Vector2 CovarianceTool::getCovariances(const SpacePoint* sp,
                                                    float zAlign, float rAlign,
                                                    float sigma) {
  Acts::Vector2 cov;
  cov[0] = ((*sp).varianceR + rAlign * rAlign) * sigma;
  cov[1] = ((*sp).varianceZ + zAlign * zAlign) * sigma;
  return cov;
}
}  // namespace Acts
