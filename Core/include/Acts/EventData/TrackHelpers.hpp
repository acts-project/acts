// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"

namespace Acts {

/// Helper function to calculate a number of track level quantities and store
/// them on the track itself
/// @note The input track needs to be mutable, so @c ReadOnly=false
/// @tparam track_container_t the track container backend
/// @tparam track_state_container_t the track state container backend
/// @tparam holder_t the holder type for the track container backends
/// @param track A mutable track proxy to operate on
template <typename track_container_t, typename track_state_container_t,
          template <typename> class holder_t>
void calculateTrackQuantities(
    Acts::TrackProxy<track_container_t, track_state_container_t, holder_t,
                     false>
        track) {
  track.chi2() = 0;
  track.nDoF() = 0;

  track.nHoles() = 0;
  track.nMeasurements() = 0;
  track.nSharedHits() = 0;
  track.nOutliers() = 0;

  for (const auto& trackState : track.trackStates()) {
    auto typeFlags = trackState.typeFlags();

    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      if (typeFlags.test(Acts::TrackStateFlag::SharedHitFlag)) {
        track.nSharedHits()++;
      }

      track.nMeasurements()++;
      track.chi2() += trackState.chi2();
      track.nDoF() += trackState.calibratedSize();
    } else if (typeFlags.test(Acts::TrackStateFlag::OutlierFlag)) {
      track.nOutliers()++;
    } else if (typeFlags.test(Acts::TrackStateFlag::HoleFlag)) {
      track.nHoles()++;
    }
  }
}
}  // namespace Acts
