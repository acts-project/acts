// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"

namespace ActsExamples {

void RefittingCalibrator::calibrate(const Acts::GeometryContext& /*gctx*/,
                                    Proxy trackState) const {
  const auto& sl =
      trackState.getUncalibratedSourceLink().get<RefittingSourceLink>();

  // Here we construct a measurement by extracting the information available
  // in the state
  Acts::visit_measurement(sl.state.calibratedSize(), [&](auto N) {
    using namespace Acts;
    constexpr int Size = decltype(N)::value;

    trackState.allocateCalibrated(Size);
    trackState.template calibrated<Size>() =
        sl.state.template calibrated<Size>();
    trackState.template calibratedCovariance<Size>() =
        sl.state.template calibratedCovariance<Size>();
  });

  trackState.setProjectorBitset(sl.state.projectorBitset());
}

}  // namespace ActsExamples
