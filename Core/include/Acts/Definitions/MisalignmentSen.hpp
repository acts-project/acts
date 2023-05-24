#pragma once

#include <iostream>
#include <unordered_map>
#include <array>

namespace Acts {
  struct Surface;
}

namespace ActsAlignment {

using SensorMisalignment = std::array<double, 6>;

namespace detail {

using namespace Acts;

struct TrackAlignmentState {
  std::unordered_map<const Surface*, SensorMisalignment> misalignments;
};

SensorMisalignment getMisalignmentParametersForSurface(const Surface* surface);

template <typename traj_t, typename parameters_t = BoundTrackParameters>
TrackAlignmentState trackAlignmentState(
    const GeometryContext& gctx,
    const Acts::MultiTrajectory<traj_t>& multiTraj,
    const size_t& entryIndex,
    const std::pair<ActsDynamicMatrix, std::unordered_map<size_t, size_t>>& globalTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces,
    const AlignmentMask& alignMask);

} 
} 

using namespace ActsAlignment::detail;
SensorMisalignment getMisalignmentParametersForSurface(const Acts::Surface* surface);
