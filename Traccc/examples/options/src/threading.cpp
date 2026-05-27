/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/threading.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <stdexcept>

namespace traccc::opts {

threading::threading() : interface("Multi-Threading Options") {

    m_desc.add_options()(
        "cpu-threads",
        boost::program_options::value(&threads)->default_value(threads),
        "The number of CPU threads to use");
}

void threading::read(const boost::program_options::variables_map &) {

    if (threads == 0) {
        throw std::invalid_argument{"Must use threads>0"};
    }
}

std::unique_ptr<configuration_printable> threading::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of CPU thread", std::to_string(threads)));

    return cat;
}

}  // namespace traccc::opts
