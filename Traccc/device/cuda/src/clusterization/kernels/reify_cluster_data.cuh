// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/edm/silicon_cluster_collection.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
namespace traccc::cuda::kernels {

/// Fill the cluster collection with the cell indices
///
/// @param disjoint_set_view The cluster/measurement index of each cell
/// @param cluster_view The collection to fill
///
__global__ void reify_cluster_data(
    vecmem::data::vector_view<const unsigned int> disjoint_set_view,
    vecmem::data::vector_view<const unsigned int> permutation_map_view,
    traccc::edm::silicon_cluster_collection::view cluster_view);

}  // namespace traccc::cuda::kernels
