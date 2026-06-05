/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/read_spacepoints.hpp"

#include "csv/read_spacepoints.hpp"
#include "read_binary.hpp"
#include "traccc/io/utils.hpp"

// System include(s).
#include <filesystem>

namespace traccc::io {

void read_spacepoints(
    edm::spacepoint_collection::host& spacepoints,
    edm::measurement_collection::host& measurements, std::size_t event,
    std::string_view directory, const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    data_format format) {

    switch (format) {
        case data_format::csv: {
            read_spacepoints(
                spacepoints, measurements,
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-hits.csv")))
                                      .native()),
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.csv")))
                                      .native()),
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurement-simhit-map.csv")))
                                      .native()),
                detector, det_desc, det_cond, format);
            break;
        }
        case data_format::binary: {
            read_spacepoints(
                spacepoints, measurements,
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-hits.dat")))
                                      .native()),
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.dat")))
                                      .native()),
                "", detector, det_desc, det_cond, format);
            break;
        }
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void read_spacepoints(
    edm::spacepoint_collection::host& spacepoints,
    edm::measurement_collection::host& measurements,
    std::string_view hit_filename, std::string_view meas_filename,
    std::string_view meas_hit_map_filename,
    const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    data_format format) {

    switch (format) {
        case data_format::csv:
            csv::read_spacepoints(spacepoints, measurements, hit_filename,
                                  meas_filename, meas_hit_map_filename,
                                  detector, det_desc, det_cond);
            break;
        case data_format::binary:
            details::read_binary_soa(spacepoints, hit_filename);
            details::read_binary_soa(measurements, meas_filename);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io
