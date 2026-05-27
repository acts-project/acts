/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/write.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/output_data.hpp"
#include "traccc/options/program_options.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdlib>

int create_binaries(const traccc::opts::detector& detector_opts,
                    const traccc::opts::input_data& input_opts,
                    const traccc::opts::output_data& output_opts,
                    std::unique_ptr<const traccc::Logger> logger) {

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Construct the detector description object.
    traccc::detector_design_description::host det_descr{host_mr};
    traccc::detector_conditions_description::host det_cond{host_mr};
    traccc::io::read_detector_description(
        det_descr, det_cond, detector_opts.detector_file,
        detector_opts.digitization_file, detector_opts.conditions_file,
        traccc::data_format::json);

    // Loop over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Read the cells from the relevant event file
        traccc::edm::silicon_cell_collection::host cells{host_mr};
        traccc::io::read_cells(cells, event, input_opts.directory,
                               logger->clone(), &det_cond, input_opts.format);

        // Write binary file
        traccc::io::write(event, output_opts.directory,
                          traccc::data_format::binary, vecmem::get_data(cells),
                          vecmem::get_data(det_descr),
                          vecmem::get_data(det_cond));

        // Read the measurements and hits from the relevant event file
        traccc::edm::measurement_collection::host measurements{host_mr};
        traccc::edm::spacepoint_collection::host spacepoints{host_mr};
        traccc::io::read_spacepoints(spacepoints, measurements, event,
                                     input_opts.directory, nullptr, nullptr,
                                     nullptr, input_opts.format);

        // Write binary file(s)
        traccc::io::write(
            event, output_opts.directory, traccc::data_format::binary,
            vecmem::get_data(spacepoints), vecmem::get_data(measurements));
    }

    return EXIT_SUCCESS;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleCreateBinaries", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::output_data output_opts;
    traccc::opts::program_options program_opts{
        "Binary File Creation",
        {detector_opts, input_opts, output_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return create_binaries(detector_opts, input_opts, output_opts,
                           logger->clone());
}
