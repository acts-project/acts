// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

template <typename index_t>
TRACCC_HOST_DEVICE inline void reduce_problem_cell(
    const edm::silicon_cell_collection::const_device& cells,
    const unsigned int cid, const unsigned int start, const unsigned int end,
    unsigned char& adjc, index_t* adjv) {
  // Some sanity check(s).
  assert(start <= end);

  // Index of the "reference cell".
  const unsigned int pos = cid + start;

  // Load the properties of the "reference cell" into local variables.
  const edm::silicon_cell reference_cell = cells.at(pos);
  const auto ref_channel0 = reference_cell.channel0();
  const auto ref_channel1 = reference_cell.channel1();
  const auto ref_module_index = reference_cell.module_index();

  /*
   * We traverse the cells backwards, starting from the current
   * cell and working back to the first, collecting adjacent cells
   * along the way.
   */
  for (unsigned int j = pos + 1; j < end; ++j) {
    const auto cell = cells.at(j);
    const auto cur_channel0 = cell.channel0();
    const auto cur_channel1 = cell.channel1();
    const auto cur_module_index = cell.module_index();

    /*
     * Since the data is sorted, we can assume that if we see a cell
     * sufficiently far away in both directions, it becomes
     * impossible for that cell to ever be adjacent to this one.
     * This is a small optimisation. It does not need abs, as channels
     * are sorted in column major.
     */
    assert((cur_channel1 >= ref_channel1) ||
           (cur_module_index != ref_module_index));
    if ((cur_channel1 > (ref_channel1 + 1)) ||
        (cur_module_index != ref_module_index)) {
      break;
    }

    /*
     * If the cell examined is adjacent to the current cell, save it
     * in the current cell's adjacency set. This uses 8-cell
     * connectivity.
     */
    if ((ref_channel0 - cur_channel0) * (ref_channel0 - cur_channel0) <= 1 &&
        (ref_channel1 - cur_channel1) * (ref_channel1 - cur_channel1) <= 1 &&
        ref_module_index == cur_module_index) {
      assert(adjc < 4);
      adjv[adjc++] = static_cast<index_t>(j - start);
    }
  }
}

}  // namespace traccc::device
