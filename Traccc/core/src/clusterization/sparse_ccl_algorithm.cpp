/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/clusterization/sparse_ccl_algorithm.hpp"

#include "traccc/clusterization/details/sparse_ccl.hpp"
#include "traccc/sanity/contiguous_on.hpp"
#include "traccc/sanity/ordered_on.hpp"
#include "traccc/utils/projections.hpp"
#include "traccc/utils/relations.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>

namespace traccc::host {

sparse_ccl_algorithm::sparse_ccl_algorithm(vecmem::memory_resource& mr,
                                           std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr(mr) {}

sparse_ccl_algorithm::output_type sparse_ccl_algorithm::operator()(
    const edm::silicon_cell_collection::const_view& cells_view,
    const detector_conditions_description::const_view& det_cond_view) const {

    // Construct the device view of the cells.
    const edm::silicon_cell_collection::const_device cells{cells_view};
    const detector_conditions_description::const_device det_cond{det_cond_view};

    // Run some sanity checks on it.
    assert(is_contiguous_on(cell_module_projection(), cells));
    assert(is_ordered_on(channel0_major_cell_order_relation(), cells));

    // Run SparseCCL to fill CCL indices.
    vecmem::vector<unsigned int> cluster_indices{cells.size(), &(m_mr.get())};
    vecmem::device_vector<unsigned int> cluster_indices_device{
        vecmem::get_data(cluster_indices)};
    const unsigned int num_clusters =
        details::sparse_ccl(cells, cluster_indices_device);

    // Create the result container.
    output_type clusters{m_mr.get()};
    clusters.resize(num_clusters);

    // Add cells to their clusters.
    for (unsigned int cell_idx = 0; cell_idx < cluster_indices.size();
         ++cell_idx) {
        clusters.cell_indices()[cluster_indices[cell_idx]].push_back(cell_idx);
    }

    // Return the clusters.
    return clusters;
}

}  // namespace traccc::host
