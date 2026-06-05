/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/magnetic_field.hpp"

// Project include(s).
#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>
#include <sstream>
#include <stdexcept>

namespace traccc::opts {

/// Helper namespace for convenience
namespace po = boost::program_options;

/// Type alias for the data format enumeration
using format_type = std::string;
/// Name of the data format option
static const char* format_option = "bfield-file-format";

magnetic_field::magnetic_field() : interface("Magnetic Field Options") {

    m_desc.add_options()("read-bfield-from-file",
                         po::bool_switch(&read_from_file),
                         "Read the magnetic field from a file");
    m_desc.add_options()("bfield-file", po::value(&file)->default_value(file),
                         "Magnetic field file");
    m_desc.add_options()(format_option,
                         po::value<format_type>()->default_value("binary"),
                         "Format of the magnetic field file");
    m_desc.add_options()("bfield-value",
                         po::value(&value)->default_value(value),
                         "Magnetic field value (when not reading from a file)");
}

void magnetic_field::read(const po::variables_map& vm) {

    // Decode the magnetic field file data format.
    if (vm.count(format_option)) {
        const std::string format_string = vm[format_option].as<format_type>();
        if (format_string == "csv") {
            format = data_format::csv;
        } else if (format_string == "binary") {
            format = data_format::binary;
        } else {
            throw std::invalid_argument("Unknown magnetic field data format");
        }
    }
}

std::unique_ptr<configuration_printable> magnetic_field::as_printable() const {

    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Read magnetic field from file", std::format("{}", read_from_file)));
    cat->add_child(
        std::make_unique<configuration_kv_pair>("Magnetic field file", file));
    std::ostringstream format_ss;
    format_ss << format;
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Magnetic field file format", format_ss.str()));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Magnetic field value", std::format("{} T", value / unit<float>::T)));

    return cat;
}

}  // namespace traccc::opts
