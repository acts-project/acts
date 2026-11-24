// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/ISurfaceMaterialAccumulator.hpp"

#include "Acts/Utilities/Diagnostics.hpp"

#include <stdexcept>

namespace Acts {

void ISurfaceMaterialAccumulator::accumulateTracks(
    State& state, const GeometryContext& gctx,
    std::span<const TrackAccumulation> tracks,
    detail::WorkerPool* /*pool*/) const {
  for (const auto& track : tracks) {
    // Support implementations that only override the deprecated method
    ACTS_PUSH_IGNORE_DEPRECATED()
    accumulate(state, gctx, *track.interactions,
               *track.surfacesWithoutAssignment);
    ACTS_POP_IGNORE_DEPRECATED()
  }
}

void ISurfaceMaterialAccumulator::accumulate(
    State& /*state*/, const GeometryContext& /*gctx*/,
    const std::vector<MaterialInteraction>& /*interactions*/,
    const std::vector<IAssignmentFinder::SurfaceAssignment>&
    /*surfacesWithoutAssignment*/) const {
  throw std::logic_error(
      "ISurfaceMaterialAccumulator implementations must override "
      "accumulateTracks");
}

}  // namespace Acts
