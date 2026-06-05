/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/read_magnetic_field.hpp"

#include "traccc/io/utils.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"

// System include(s).
#include <format>
#include <fstream>
#include <stdexcept>

namespace traccc::io {
namespace binary {

void read_magnetic_field(magnetic_field& bfield, std::string_view filename,
                         std::unique_ptr<const Logger> ilogger) {

    // Set up a local logger.
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Open the file.
    std::ifstream ifile{get_absolute_path(filename),
                        std::ios::binary | std::ios::in};
    if (!ifile.is_open()) {
        throw std::invalid_argument(
            std::format("Failed to open magnetic field file: {}", filename));
    }

    // Construct/fill the magnetic field from the file.
    TRACCC_INFO("Reading magnetic field from file: " << filename);
    bfield.set(covfie::field<host::inhom_bfield_backend_t<traccc::scalar>>(
        covfie::field<host::inhom_io_bfield_backend_t<traccc::scalar>>(ifile)));
}

}  // namespace binary

void read_magnetic_field(magnetic_field& bfield, std::string_view filename,
                         data_format format,
                         std::unique_ptr<const Logger> ilogger) {

    switch (format) {
        case data_format::binary:
            binary::read_magnetic_field(bfield, filename, std::move(ilogger));
            break;
        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io
