// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
