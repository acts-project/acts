// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../../utils/thread_id.hpp"
#include "reify_cluster_data.cuh"

// Project include(s).
#include "traccc/clusterization/device/reify_cluster_data.hpp"

namespace traccc::cuda::kernels {
__global__ void reify_cluster_data(
    vecmem::data::vector_view<const unsigned int> disjoint_set_view,
    vecmem::data::vector_view<const unsigned int> permutation_map_view,
    traccc::edm::silicon_cluster_collection::view cluster_view) {
  device::reify_cluster_data(details::thread_id1{}.getGlobalThreadId(),
                             disjoint_set_view, permutation_map_view,
                             cluster_view);
}
}  // namespace traccc::cuda::kernels
