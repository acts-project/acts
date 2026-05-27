/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/clusterization/measurement_sorting_algorithm.hpp"

// System include(s).
#include <algorithm>
#include <numeric>

namespace traccc::host {

measurement_sorting_algorithm::measurement_sorting_algorithm(
    vecmem::memory_resource& mr, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr(mr) {}

measurement_sorting_algorithm::output_type
measurement_sorting_algorithm::operator()(
    const edm::measurement_collection::const_view& measurements_view) const {

    // Create a device container on top of the view.
    const edm::measurement_collection::const_device measurements{
        measurements_view};

    // Create a vector of measurement indices, which would be sorted.
    vecmem::vector<unsigned int> indices(measurements.size(), &(m_mr.get()));
    std::iota(indices.begin(), indices.end(), 0u);

    // Sort the indices according to the measurements.
    std::sort(indices.begin(), indices.end(),
              [&](unsigned int lhs, unsigned int rhs) {
                  return measurements.at(lhs) < measurements.at(rhs);
              });

    // Fill an output container with the sorted measurements.
    edm::measurement_collection::host result{m_mr.get()};
    for (unsigned int i : indices) {
        result.push_back(measurements.at(i));
    }

    // Return the sorted measurements.
    return result;
}

}  // namespace traccc::host
