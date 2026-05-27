/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/program_options.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <cstdlib>
#include <iostream>

namespace traccc::opts {

program_options::program_options(
    std::string_view description,
    const std::vector<std::reference_wrapper<interface> >& options, int argc,
    char* argv[], std::unique_ptr<const traccc::Logger> ilogger)
    : m_desc(std::string{description}) {

    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Add all of the option groups.
    for (const interface& opt : options) {
        m_desc.add(opt.options());
    }
    // Add a help option.
    m_desc.add_options()("help,h", "Print this help message");

    // Parse the command line options.
    boost::program_options::variables_map vm;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, m_desc), vm);

    // Print a help message if the user asked for it.
    if (vm.count("help")) {
        std::cout << m_desc << std::endl;
        std::exit(0);
    }

    // Handle any and all errors.
    try {
        boost::program_options::notify(vm);
    } catch (const std::exception& ex) {
        TRACCC_FATAL("Couldn't interpret command line options because of: "
                     << ex.what() << "; " << m_desc);
        std::exit(1);
    }

    // Read / post-process the options.
    for (interface& opt : options) {
        opt.read(vm);
    }

    // Tell the user what's happening.
    TRACCC_INFO("\nRunning " << description);

    configuration_list cl;

    for (const auto& opt : options) {
        cl.add_child(opt.get().as_printable());
    }

    TRACCC_INFO("\n" << cl.print() << "\n");
}

}  // namespace traccc::opts
