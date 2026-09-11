/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::details {

/// Implemementation of SparseCCL, following
/// [DOI: 10.1109/DASIP48288.2019.9049184]
///
/// Requires cells to be sorted in column major

/// Find root of the tree for entry @param e
///
/// @param labels an equivalance table
///
/// @return the root of @param e
///
TRACCC_HOST_DEVICE inline unsigned int find_root(
    const vecmem::device_vector<unsigned int>& labels, unsigned int e);

/// Create a union of two entries @param e1 and @param e2
///
/// @param labels an equivalance table
///
/// @return the rleast common ancestor of the entries
///
TRACCC_HOST_DEVICE inline unsigned int make_union(
    vecmem::device_vector<unsigned int>& labels, unsigned int e1,
    unsigned int e2);

/// Helper method to find adjacent cells
///
/// @param a the first cell
/// @param b the second cell
///
/// @return boolan to indicate 8-cell connectivity
///
template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline bool is_adjacent(const edm::silicon_cell<T1>& a,
                                           const edm::silicon_cell<T2>& b);

/// Helper method to find define distance,
/// does not need abs, as channels are sorted in
/// column major
///
/// @param a the first cell
/// @param b the second cell
///
/// @return boolan to indicate !8-cell connectivity
///
template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline bool is_far_enough(const edm::silicon_cell<T1>& a,
                                             const edm::silicon_cell<T2>& b);

/// Sparce CCL algorithm
///
/// @param cells is the cell collection
/// @param labels is the vector of the output indices (to which cluster a cell
///               belongs to)
/// @return number of clusters
///
TRACCC_HOST_DEVICE inline unsigned int sparse_ccl(
    const edm::silicon_cell_collection::const_device& cells,
    vecmem::device_vector<unsigned int>& labels);

}  // namespace traccc::details

// Include the implementation.
#include "traccc/clusterization/impl/sparse_ccl.ipp"
