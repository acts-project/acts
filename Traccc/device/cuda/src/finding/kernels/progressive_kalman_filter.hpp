// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/finding/device/progressive_kalman_filter.hpp"
#include "traccc/finding/finding_config.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

template <typename propagator_t>
void progressive_kalman_filter(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const typename propagator_t::stepper_type::magnetic_field_type& field_data,
    const vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>& surfaces_view,
    const device::progressive_kalman_filter_payload& payload);

}  // namespace traccc::cuda
