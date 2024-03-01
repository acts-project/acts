// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

namespace Acts {
class ConstVectorMultiTrajectory;
class VectorMultiTrajectory;
}  // namespace Acts

namespace ActsExamples {

struct RefittingCalibrator {
  using Proxy = Acts::VectorMultiTrajectory::TrackStateProxy;
  using ConstProxy = Acts::ConstVectorMultiTrajectory::ConstTrackStateProxy;

  struct RefittingSourceLink {
    ConstProxy state;

    Acts::GeometryIdentifier geometryId() const {
      return state.referenceSurface().geometryId();
    }
  };

  void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink, Proxy trackState) const;
};

}  // namespace ActsExamples
