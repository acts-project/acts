/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/read_measurements.hpp"

#include "csv/read_measurements.hpp"
#include "read_binary.hpp"
#include "traccc/io/utils.hpp"

// System include(s).
#include <filesystem>

namespace traccc::io {

std::vector<measurement_id_type> read_measurements(
    edm::measurement_collection::host& measurements, std::size_t event,
    std::string_view directory, const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    const bool sort_measurements, data_format format) {

    switch (format) {
        case data_format::csv: {
            return read_measurements(
                measurements,
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.csv")))
                                      .native()),
                detector, det_desc, det_cond, sort_measurements, format);
        }
        case data_format::binary: {
            return read_measurements(
                measurements,
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.dat")))
                                      .native()),
                detector, det_desc, det_cond, sort_measurements, format);
        }
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

std::vector<measurement_id_type> read_measurements(
    edm::measurement_collection::host& measurements, std::string_view filename,
    const traccc::host_detector* detector,
    const traccc::detector_design_description::host* det_desc,
    const traccc::detector_conditions_description::host* det_cond,
    const bool sort_measurements, data_format format) {

    switch (format) {
        case data_format::csv:
            return csv::read_measurements(measurements, filename, detector,
                                          det_desc, det_cond,
                                          sort_measurements);
        case data_format::binary:
            details::read_binary_soa(measurements, filename);
            return {};
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io
