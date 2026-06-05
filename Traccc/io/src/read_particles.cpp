/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/read_particles.hpp"

#include "csv/read_particles.hpp"
#include "traccc/io/utils.hpp"

// System include(s).
#include <filesystem>

namespace traccc::io {

void read_particles(particle_collection_types::host& particles,
                    std::size_t event, std::string_view directory,
                    data_format format, std::string_view filename_postfix) {

    switch (format) {
        case data_format::csv:
            read_particles(
                particles,
                get_absolute_path(
                    (std::filesystem::path(directory) /
                     std::filesystem::path(get_event_filename(
                         event, std::string{filename_postfix} + ".csv")))
                        .native()),
                format);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void read_particles(particle_collection_types::host& particles,
                    std::string_view filename, data_format format) {

    switch (format) {
        case data_format::csv:
            csv::read_particles(particles, filename);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void read_particles(particle_container_types::host& particles,
                    edm::measurement_collection::host& measurements,
                    std::size_t event, std::string_view directory,
                    const traccc::host_detector* detector, data_format format,
                    std::string_view filename_postfix) {

    switch (format) {
        case data_format::csv:
            read_particles(
                particles, measurements,
                get_absolute_path(
                    (std::filesystem::path(directory) /
                     std::filesystem::path(get_event_filename(
                         event, std::string{filename_postfix} + ".csv")))
                        .native()),
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
                detector, format);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void read_particles(particle_container_types::host& particles,
                    edm::measurement_collection::host& measurements,
                    std::string_view particles_file, std::string_view hits_file,
                    std::string_view measurements_file,
                    std::string_view hit_map_file,
                    const traccc::host_detector* detector, data_format format) {

    switch (format) {
        case data_format::csv:
            csv::read_particles(particles, measurements, particles_file,
                                hits_file, measurements_file, hit_map_file,
                                detector);
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io
