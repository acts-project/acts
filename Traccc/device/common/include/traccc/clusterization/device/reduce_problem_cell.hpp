// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"

namespace traccc::device {

/// Function for looking for adjacent cells. The cell ids will range from 0 to
/// max_cells_per_partition and the number of blocks will equal the number of
/// partitions, hence checking all cells.
///
/// @param[in] cells    Collection of cells
/// @param[in] cid      Current cell id
/// @param[in] start    Current partition start point
/// @param[in] end      Current partition end point
/// @param[out] ajc     Number of adjacent cells
/// @param[out] ajv     Indices of adjacent cells
///
template <typename index_t>
TRACCC_HOST_DEVICE inline void reduce_problem_cell(
    const edm::silicon_cell_collection::const_device& cells, unsigned int cid,
    unsigned int start, unsigned int end, unsigned char& adjc, index_t* adjv);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/clusterization/device/impl/reduce_problem_cell.ipp"
