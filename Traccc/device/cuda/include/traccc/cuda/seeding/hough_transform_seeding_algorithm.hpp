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
#include "traccc/seeding/device/hough_transform_seeding_algorithm.hpp"

namespace traccc::cuda {

/// Main algorithm for performing the track seeding on an NVIDIA GPU
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class hough_transform_seeding_algorithm	: public device::hough_transform_seeding_algorithm,
                                  	  public cuda::algorithm_base {
 public:
  /// Constructor for the seed finding algorithm
  ///
  /// @param finder_config The seed finding configuration
  /// @param grid_config   The spacepoint grid forming configuration
  /// @param filter_config The seed filtering configuration
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param str The CUDA stream to perform the operations in
  /// @param logger The logger instance to use for messaging
  /// @param await_func The function to use for synchronizing events
  ///
  hough_transform_seeding_algorithm(const memory_resource& mr,
				    const vecmem::copy& copy,
				    const stream_wrapper& str,
				    const hough_seeding_config& seeding_config,
				    std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
				    await_function_type await_func = await_sync_event);

  private:
  void define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const;
  void precompute_spacepoint_variables_kernel(const precompute_spacepoint_variables_payload& payload) const;
  void fill_hough_plane_kernel(const fill_hough_planes_kernel_payload& payload) const;
  void count_maxima_kernel(const count_maxima_payload& payload) const;
  void find_maxima_kernel(const find_maxima_payload& payload) const;
  void take_picture_kernel(const take_picture_payload& payload) const;
  void evaluate_nn_kernel(const evaluate_nn_payload& payload) const;
  void count_contributing_hits_kernel(const count_contributing_hits_payload& payload) const;

  /// Wait for outstanding work on the algorithm stream or queue.
  void synchronize() const override;

};

}  // namespace traccc::cuda
