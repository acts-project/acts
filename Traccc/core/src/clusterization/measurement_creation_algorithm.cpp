/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/clusterization/measurement_creation_algorithm.hpp"

#include "traccc/clusterization/details/measurement_creation.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <cassert>

namespace traccc::host {

measurement_creation_algorithm::measurement_creation_algorithm(
    vecmem::memory_resource &mr, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr(mr) {}

measurement_creation_algorithm::output_type
measurement_creation_algorithm::operator()(
    const edm::silicon_cell_collection::const_view &cells_view,
    const edm::silicon_cluster_collection::const_view &clusters_view,
    const detector_design_description::const_view &dmd_view,
    const detector_conditions_description::const_view &dcd_view) const {

    // Create device containers for the input variables.
    const edm::silicon_cell_collection::const_device cells{cells_view};
    const edm::silicon_cluster_collection::const_device clusters{clusters_view};
    const detector_design_description::const_device det_descr{dmd_view};
    const detector_conditions_description::const_device det_cond{dcd_view};

    // Create the result object.
    output_type result(m_mr.get());
    result.resize(clusters.size());
    edm::measurement_collection::device measurements{vecmem::get_data(result)};

    // Process the clusters one-by-one.
    for (decltype(clusters)::size_type i = 0; i < clusters.size(); ++i) {

        // Get the cluster and measurement.
        const edm::silicon_cluster cluster = clusters.at(i);
        edm::measurement measurement = measurements.at(i);

        // A security check.
        assert(cluster.cell_indices().empty() == false);

        // Fill measurement from cluster
        details::fill_measurement(measurement, cluster, i, cells, det_descr,
                                  det_cond);
    }

    return result;
}

}  // namespace traccc::host
