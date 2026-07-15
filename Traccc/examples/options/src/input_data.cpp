/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/input_data.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>
#include <sstream>
#include <stdexcept>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

/// Type alias for the data format enumeration
using data_format_type = std::string;
/// Name of the data format option
static const char *data_format_option = "input-data-format";

input_data::input_data() : interface("Input Data Options") {

    m_desc.add_options()(
        "use-acts-geom-source",
        po::value(&use_acts_geom_source)->default_value(use_acts_geom_source),
        "Use acts geometry source");
    m_desc.add_options()(data_format_option,
                         po::value<data_format_type>()->default_value("csv"),
                         "Format of the input file(s)");
    m_desc.add_options()("input-directory",
                         po::value(&directory)->default_value(directory),
                         "Directory holding the input files");
    m_desc.add_options()("input-events",
                         po::value(&events)->default_value(events),
                         "Number of input events to process");
    m_desc.add_options()("input-skip", po::value(&skip)->default_value(skip),
                         "Number of input events to skip");
}

void input_data::read(const po::variables_map &vm) {

    // Decode the input data format.
    if (vm.count(data_format_option)) {
        const std::string input_format_string =
            vm[data_format_option].as<data_format_type>();
        if (input_format_string == "csv") {
            format = data_format::csv;
        } else if (input_format_string == "binary") {
            format = data_format::binary;
        } else if (input_format_string == "json") {
            format = data_format::json;
        } else {
            throw std::invalid_argument("Unknown input data format");
        }
    }
}

std::unique_ptr<configuration_printable> input_data::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Use ACTS geometry source", std::format("{}", use_acts_geom_source)));
    std::ostringstream format_ss;
    format_ss << format;
    cat->add_child(std::make_unique<configuration_kv_pair>("Input data format",
                                                           format_ss.str()));
    cat->add_child(
        std::make_unique<configuration_kv_pair>("Input directory", directory));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of input events", std::to_string(events)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of skipped events", std::to_string(skip)));

    return cat;
}

}  // namespace traccc::opts
