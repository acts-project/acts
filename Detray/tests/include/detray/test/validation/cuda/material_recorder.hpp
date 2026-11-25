// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/core/detector.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/validation/material_recording.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// System include
#include <string_view>

namespace detray::validation::cuda {

/// Launch the material validation kernel
///
/// @param[in] det_view the detector vecmem view
/// @param[in] cfg the propagation configuration
/// @param[in] tracks_view the initial track parameter of every test track
/// @param[out] track_mat_view the accumulated material per track
template <concepts::detector detector_t>
void record_material(
    typename detector_t::view_type det_view, const propagation::config &cfg,
    vecmem::data::vector_view<
        free_track_parameters<typename detector_t::algebra_type>> &tracks_view,
    vecmem::data::vector_view<track_material> &track_mat_view,
    vecmem::data::jagged_vector_view<material_record> &mat_steps_view);

/// Prepare data for device material recording run
struct material_recorder {
  /// Test name suffix
  static constexpr std::string_view name = "cuda";

  /// Run the propagation that records the material
  /// @param host_mr memory resource for host allocations
  /// @param dev_mr memory resource for device allocations
  /// @param det the detector for which the material should be recorded
  /// @param cfg propagation configuration
  /// @param tracks input tracks for which to record the material
  /// @param capacities capacities for the data recording buffer
  /// @returns the accumulated track material and material steps per track
  template <concepts::detector detector_t>
  auto operator()(
      vecmem::memory_resource *host_mr, vecmem::memory_resource *dev_mr,
      const detector_t &det, const propagation::config &cfg,
      const vecmem::vector<
          free_track_parameters<typename detector_t::algebra_type>> &tracks,
      const std::vector<std::size_t> &capacities) const {
    using algebra_t = typename detector_t::algebra_type;
    using track_t = free_track_parameters<algebra_t>;
    using material_record_t = material_record;
    using track_material_t = track_material;

    // Helper object for performing memory copies (to CUDA devices)
    vecmem::cuda::copy cuda_cpy;

    // Copy the detector to device and get its view
    const auto det_buffer = detray::get_buffer(det, *dev_mr, cuda_cpy);
    auto det_view = detray::get_data(det_buffer);

    // Move the track parameters to device
    auto tracks_buffer = cuda_cpy.to(vecmem::get_data(tracks), *dev_mr,
                                     vecmem::copy::type::host_to_device);
    vecmem::data::vector_view<track_t> tracks_view =
        vecmem::get_data(tracks_buffer);

    vecmem::data::vector_buffer<track_material_t> track_mat_buffer(
        static_cast<unsigned int>(tracks.size()), *dev_mr,
        vecmem::data::buffer_type::fixed_size);
    cuda_cpy.setup(track_mat_buffer)->wait();
    auto track_mat_view = vecmem::get_data(track_mat_buffer);

    // Buffer for the material parameters at every surface per track
    vecmem::data::jagged_vector_buffer<material_record_t> mat_steps_buffer(
        capacities, *dev_mr, host_mr, vecmem::data::buffer_type::resizable);
    cuda_cpy.setup(mat_steps_buffer)->wait();
    auto mat_steps_view = vecmem::get_data(mat_steps_buffer);

    // Run the material tracing on device
    detray::validation::cuda::record_material<detector_t>(
        det_view, cfg, tracks_view, track_mat_view, mat_steps_view);

    // Get the results back to the host and pass them on to be checked
    vecmem::vector<track_material_t> track_mat_vec(host_mr);
    cuda_cpy(track_mat_buffer, track_mat_vec)->wait();

    vecmem::jagged_vector<material_record_t> mat_steps(host_mr);
    cuda_cpy(mat_steps_buffer, mat_steps)->wait();

    return std::make_tuple(track_mat_vec, mat_steps);
  }
};

}  // namespace detray::validation::cuda
