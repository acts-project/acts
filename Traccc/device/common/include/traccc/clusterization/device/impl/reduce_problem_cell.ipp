/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/clusterization/details/sparse_ccl.hpp"
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

  // Load the "reference cell" into a local variable.
  const edm::silicon_cell reference_cell = cells.at(pos);

  /*
   * We traverse the cells backwards, starting from the current
   * cell and working back to the first, collecting adjacent cells
   * along the way.
   */
  for (unsigned int j = pos + 1; j < end; ++j) {
    /*
     * Since the data is sorted, we can assume that if we see a cell
     * sufficiently far away in both directions, it becomes
     * impossible for that cell to ever be adjacent to this one.
     * This is a small optimisation.
     */
    if (traccc::details::is_far_enough(cells.at(j), reference_cell)) {
      break;
    }

    /*
     * If the cell examined is adjacent to the current cell, save it
     * in the current cell's adjacency set.
     */
    if (traccc::details::is_adjacent(reference_cell, cells.at(j))) {
      assert(adjc < 4);
      adjv[adjc++] = static_cast<index_t>(j - start);
    }
  }
}

}  // namespace traccc::device
