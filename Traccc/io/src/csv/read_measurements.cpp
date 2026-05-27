/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_measurements.hpp"

#include "traccc/io/csv/make_measurement_edm.hpp"
#include "traccc/io/csv/make_measurement_reader.hpp"

// System include(s).
#include <numeric>
#include <ranges>

namespace traccc::io::csv {

std::vector<measurement_id_type> read_measurements(
    edm::measurement_collection::host& measurements, std::string_view filename,
    const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    const bool do_sort) {

    // Construct the measurement reader object.
    auto reader = make_measurement_reader(filename);

    // For Acts data, build a map of acts->detray geometry IDs
    std::map<geometry_id, geometry_id> acts_to_detray_id;
    std::map<geometry_id, std::size_t>
        geometry_id_to_detector_description_index;

    if (detector) {
        host_detector_visitor<detector_type_list>(
            *detector, [&]<typename detector_t>(const detector_t::host& det) {
                for (const auto& surface_desc : det.surfaces()) {
                    acts_to_detray_id[surface_desc.source] =
                        surface_desc.identifier().value();
                }
            });
    }

    if (det_cond) {

        for (std::size_t i = 0; i < det_cond->geometry_id().size(); ++i) {
            geometry_id_to_detector_description_index
                [det_cond->geometry_id().at(i).value()] =
                    det_cond->module_to_design_id().at(i);
        }
    }

    // Read the measurements from the input file.
    csv::measurement iomeas;
    while (reader.read(iomeas)) {

        // Construct the measurement object.
        measurements.resize(measurements.size() + 1u);
        edm::measurement meas = measurements.at(measurements.size() - 1u);
        make_measurement_edm(
            iomeas, meas, (detector == nullptr ? nullptr : &acts_to_detray_id),
            det_desc,
            (det_cond == nullptr ? nullptr
                                 : &geometry_id_to_detector_description_index));
    }

    // Contains the index of the new position at the entry of the old position
    std::vector<measurement_id_type> new_idx_map(measurements.size());
    if (do_sort) {
        // Remeber index locations
        std::vector<unsigned int> idx(measurements.size());
        std::iota(idx.begin(), idx.end(), 0u);

        // Sort the indices the way the measurements will be sorted
        // https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
        std::ranges::sort(idx.begin(), idx.end(),
                          [&measurements](unsigned int i, unsigned int j) {
                              return measurements[i] < measurements[j];
                          });

        // Create a sorted measurement collection.
        edm::measurement_collection::host sorted_measurements(
            measurements.resource());
        sorted_measurements.resize(measurements.size());

        // Map the indices to the new positions. While creating a sorted
        // measurement collection.
        for (std::size_t i = 0u; i < idx.size(); ++i) {
            new_idx_map[idx[i]] = static_cast<measurement_id_type>(i);
            sorted_measurements.at(i) = measurements.at(idx[i]);
        }

        // Override the measurements with the sorted ones.
        measurements = sorted_measurements;
    }

    return new_idx_map;
}

}  // namespace traccc::io::csv
