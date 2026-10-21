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

namespace {
/// The accumulator whose default implementation is currently forwarding to
/// the other one, to catch implementations that override neither of them
thread_local const ISurfaceMaterialAccumulator* tlForwarding = nullptr;

class ForwardingGuard {
 public:
  explicit ForwardingGuard(const ISurfaceMaterialAccumulator* accumulator)
      : m_previous(tlForwarding) {
    if (tlForwarding == accumulator) {
      throw std::logic_error(
          "ISurfaceMaterialAccumulator implementations must override "
          "accumulateTracks");
    }
    tlForwarding = accumulator;
  }
  ForwardingGuard(const ForwardingGuard&) = delete;
  ForwardingGuard& operator=(const ForwardingGuard&) = delete;
  ~ForwardingGuard() { tlForwarding = m_previous; }

 private:
  const ISurfaceMaterialAccumulator* m_previous;
};
}  // namespace

void ISurfaceMaterialAccumulator::accumulateTracks(
    State& state, const GeometryContext& gctx,
    std::span<const TrackAccumulation> tracks,
    detail::WorkerPool* /*pool*/) const {
  ForwardingGuard guard(this);
  for (const auto& track : tracks) {
    // Support implementations that only override the deprecated method
    ACTS_PUSH_IGNORE_DEPRECATED()
    accumulate(state, gctx, *track.interactions,
               *track.surfacesWithoutAssignment);
    ACTS_POP_IGNORE_DEPRECATED()
  }
}

void ISurfaceMaterialAccumulator::accumulate(
    State& state, const GeometryContext& gctx,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>&
        surfacesWithoutAssignment) const {
  ForwardingGuard guard(this);
  TrackAccumulation track{&interactions, &surfacesWithoutAssignment};
  accumulateTracks(state, gctx, std::span(&track, 1), nullptr);
}

}  // namespace Acts
