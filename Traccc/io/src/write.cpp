/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/write.hpp"

#include "csv/write_cells.hpp"
#include "json/write_digitization_config.hpp"
#include "obj/write_seeds.hpp"
#include "obj/write_spacepoints.hpp"
#include "obj/write_tracks.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/utils.hpp"
#include "write_binary.hpp"

// System include(s).
#include <filesystem>
#include <stdexcept>

namespace traccc::io {

void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           traccc::edm::silicon_cell_collection::const_view cells,
           traccc::detector_design_description::const_view det_desc_view,
           traccc::detector_conditions_description::const_view det_cond_view,
           bool use_acts_geometry_id) {

    switch (format) {
        case data_format::binary:
            details::write_binary_soa(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-cells.dat")))
                                      .native()),
                traccc::edm::silicon_cell_collection::const_device{cells});
            break;
        case data_format::csv:
            csv::write_cells(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-cells.csv")))
                                      .native()),
                cells, det_desc_view, det_cond_view, use_acts_geometry_id);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::spacepoint_collection::const_view spacepoints,
           edm::measurement_collection::const_view measurements) {

    switch (format) {
        case data_format::binary:
            details::write_binary_soa(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-hits.dat")))
                                      .native()),
                edm::spacepoint_collection::const_device{spacepoints});
            details::write_binary_soa(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.dat")))
                                      .native()),
                edm::measurement_collection::const_device{measurements});
            break;
        case data_format::obj:
            obj::write_spacepoints(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-spacepoints.obj")))
                                      .native()),
                spacepoints);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::measurement_collection::const_view measurements) {

    switch (format) {
        case data_format::binary:
            details::write_binary_soa(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-measurements.dat")))
                                      .native()),
                edm::measurement_collection::const_device{measurements});
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void write(std::size_t event, std::string_view directory,
           traccc::data_format format, edm::seed_collection::const_view seeds,
           edm::spacepoint_collection::const_view spacepoints) {

    switch (format) {
        case data_format::obj:
            obj::write_seeds(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(
                                       get_event_filename(event, "-seeds.obj")))
                                      .native()),
                seeds, spacepoints);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void write(std::size_t event, std::string_view directory,
           traccc::data_format format,
           edm::track_container<default_algebra>::const_view tracks,
           const traccc::host_detector& detector) {

    switch (format) {
        case data_format::obj:
            obj::write_tracks(
                get_absolute_path((std::filesystem::path(directory) /
                                   std::filesystem::path(get_event_filename(
                                       event, "-track-candidates.obj")))
                                      .native()),
                tracks, detector);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void write(std::string_view filename, data_format format,
           const digitization_config& config) {

    switch (format) {
        case data_format::json:
            json::write_digitization_config(get_absolute_path(filename),
                                            config);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io
