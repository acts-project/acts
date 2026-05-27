/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_spacepoints.hpp"

#include "read_measurements.hpp"
#include "traccc/io/csv/make_hit_reader.hpp"
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <ranges>
#include <stdexcept>
#include <unordered_map>

namespace traccc::io::csv {

void read_spacepoints(
    edm::spacepoint_collection::host& spacepoints,
    edm::measurement_collection::host& measurements,
    std::string_view hit_filename, std::string_view meas_filename,
    std::string_view meas_hit_map_filename,
    const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    const bool sort_measurements) {

    // Read all measurements.
    const std::vector<measurement_id_type> new_idx_map =
        read_measurements(measurements, meas_filename, detector, det_desc,
                          det_cond, sort_measurements);

    // Measurement hit id reader
    auto mhid_reader =
        io::csv::make_measurement_hit_id_reader(meas_hit_map_filename);
    std::unordered_map<std::size_t, traccc::io::csv::measurement_hit_id>
        measurement_hit_ids;
    traccc::io::csv::measurement_hit_id io_mh_id;
    while (mhid_reader.read(io_mh_id)) {
        if (sort_measurements) {
            io_mh_id.measurement_id = new_idx_map[io_mh_id.measurement_id];
        }
        measurement_hit_ids.insert({io_mh_id.hit_id, io_mh_id});
    }

    // Construct the hit reader object.
    auto hit_reader = make_hit_reader(hit_filename);

    // Read the hits from the input file.
    hit iohit;
    while (hit_reader.read(iohit)) {

        // Find the index of the measurement that this hit/spacepoint belongs
        // to. Which may not be valid, as some simulated hits are not associated
        // with a measurement.
        const auto measurement_id_it =
            measurement_hit_ids.find(spacepoints.size());
        const unsigned int measurement_index =
            (measurement_id_it != measurement_hit_ids.end())
                ? static_cast<unsigned int>(
                      measurement_id_it->second.measurement_id)
                : static_cast<unsigned int>(-1);

        // Create a new spacepoint for the SoA container.
        spacepoints.push_back(
            {measurement_index,
             edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
             {iohit.tx, iohit.ty, iohit.tz},
             0.f,
             0.f});
    }
}

}  // namespace traccc::io::csv
