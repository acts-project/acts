/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "kalman_fitting_telescope_test.hpp"

namespace traccc {

/// Combinatorial Kalman Finding Test with Sparse tracks
class CkfSparseTrackTelescopeTests : public KalmanFittingTelescopeTests {};

/// Combinatorial Kalman Finding Test with Identical tracks
class CkfCombinatoricsTelescopeTests : public KalmanFittingTelescopeTests {};

/// Combinatorial Kalman Finding Test with Identical tracks (CPU)
class CpuCkfCombinatoricsTelescopeTests
    : public CkfCombinatoricsTelescopeTests {};

/// Combinatorial Kalman Finding Test with Identical tracks (CUDA)
class CudaCkfCombinatoricsTelescopeTests
    : public CkfCombinatoricsTelescopeTests {};

}  // namespace traccc
