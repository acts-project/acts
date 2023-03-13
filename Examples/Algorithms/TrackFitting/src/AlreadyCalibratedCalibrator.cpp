// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/AlreadyCalibratedCalibrator.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"

namespace ActsExamples {

void AlreadyCalibratedCalibrator::calibrate(
    const Acts::GeometryContext& /*gctx*/, Proxy trackState) const {
  const ActsExamples::IndexSourceLink& sourceLink =
      trackState.getUncalibratedSourceLink()
          .get<ActsExamples::IndexSourceLink>();

  ConstProxy old = callibratedStates.at(sourceLink.index());

  // Here we construct a measurement by extracting the information available
  // in the state
  Acts::visit_measurement(old.calibratedSize(), [&](auto N) {
    using namespace Acts;
    constexpr int Size = decltype(N)::value;

    std::array<Acts::BoundIndices, Size> indices;
    auto it = indices.begin();

    for (auto eBound : {eBoundLoc0, eBoundLoc1, eBoundPhi, eBoundTheta,
                        eBoundQOverP, eBoundTime}) {
      bool haveIndex = old.effectiveProjector().col(eBound).sum() == 0;
      if (haveIndex) {
        *it = eBound;
        ++it;
      }
    }

    Acts::Measurement<BoundIndices, Size> measurement(
        SourceLink(sourceLink), indices, old.template calibrated<Size>(),
        old.template calibratedCovariance<Size>());

    trackState.allocateCalibrated(Size);
    trackState.setCalibrated<Size>(measurement);
  });
}
};  // namespace ActsExamples
