// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Traccc/TrackConversion.hpp"

namespace ActsPlugins::TracccBoundTrackParameters {

Acts::BoundVector convertBoundVector(const BoundTrackParameters& params) {
  Acts::BoundVector vector;
  for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
    vector[i] = traccc::getter::element(params.vector(), i, 0u);
  }
  return vector;
}

Acts::BoundMatrix convertBoundCovariance(const BoundTrackParameters& params) {
  Acts::BoundMatrix covariance;
  for (unsigned int i = 0; i < Acts::eBoundSize; ++i) {
    for (unsigned int j = 0; j < Acts::eBoundSize; ++j) {
      covariance(i, j) = traccc::getter::element(params.covariance(), i, j);
    }
  }
  return covariance;
}

bool isDegenerate(const BoundTrackParameters& params) {
  return params.bound_local()[0] == 0 && params.bound_local()[1] == 0 &&
         params.phi() == 0 && params.theta() == 0 && params.qop() == 0 &&
         params.time() == 0;
}

std::ostream& operator<<(std::ostream& os,
                         const TrackConversionStatistics& stats) {
  os << "converted " << stats.nConverted << " / " << stats.nInput
     << " tracks (dropped: " << stats.nFailedFit << " failed fit, "
     << stats.nMissingSurface << " missing surface, " << stats.nDegenerate
     << " degenerate parameters)";
  return os;
}

}  // namespace ActsPlugins::Traccc
