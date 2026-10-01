// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/finding/device/find_tracks.hpp"
#include "traccc/finding/finding_config.hpp"

namespace traccc::cuda {

template <detray::concepts::detector detector_t>
void find_tracks(const dim3& grid_size, const dim3& block_size,
                 std::size_t shared_mem_size, const cudaStream_t& stream,
                 const finding_config& cfg,
                 const typename detector_t::const_view_type& det,
                 const device::find_tracks_payload& payload);
}  // namespace traccc::cuda
