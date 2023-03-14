// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace ActsExamples {

struct RefittingCalibrator {
  using Proxy =
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy;
  using ConstProxy = Acts::MultiTrajectory<
      Acts::ConstVectorMultiTrajectory>::ConstTrackStateProxy;

  struct RefittingSourceLink {
    ConstProxy state;

    Acts::GeometryIdentifier geometryId() const {
      return state.referenceSurface().geometryId();
    }
  };

  void calibrate(const Acts::GeometryContext& /*gctx*/, Proxy trackState) const;
};

}  // namespace ActsExamples
