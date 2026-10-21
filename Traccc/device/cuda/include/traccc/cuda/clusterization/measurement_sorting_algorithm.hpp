// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"
#include "traccc/cuda/utils/await.hpp"

// Project include(s).
#include "traccc/clusterization/device/measurement_sorting_algorithm.hpp"

namespace traccc::cuda {

/// Algorithm sorting the reconstructed measurements
///
/// The track finding algorithm expects measurements belonging to a single
/// detector module to be consecutive in memory. But
/// @c traccc::cuda::clusterization_algorithm does not (currently) produce the
/// measurements in such an ordered state. This is where this algorithm comes
/// to the rescue.
///
class measurement_sorting_algorithm
    : public device::measurement_sorting_algorithm,
      public cuda::algorithm_base {
 public:
  /// Constructor for the algorithm
  ///
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use in the algorithm
  /// @param str The CUDA stream to schedule the measurement sorting in
  /// @param logger The logger to use in the algorithm
  /// @param await_func The function to use for synchronizing events
  ///
  measurement_sorting_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      const stream_wrapper& str,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
      await_function_type await_func = await_sync_event);

 private:
  /// Run the measurement-sorting kernels on the CUDA stream.
  void sorting_kernel(
      const measurement_sorting_kernel_payload& payload) const override;

  /// Wait for outstanding work on the algorithm stream or queue.
  void synchronize() const override;
};  // class measurement_sorting_algorithm

}  // namespace traccc::cuda
