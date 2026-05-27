/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// core
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/utils/propagation.hpp"

// algorithms
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/sycl/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/sycl/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/sycl/utils/make_magnetic_field.hpp"

// io
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/utils.hpp"

// performance
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/soa_comparator.hpp"
#include "traccc/performance/timer.hpp"

// options
#include "traccc/options/accelerator.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_seeding.hpp"

// examples
#include "traccc/examples/make_magnetic_field.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

int seq_run(const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Creating sycl queue object
    vecmem::sycl::queue_wrapper vecmem_queue;
    traccc::sycl::queue_wrapper traccc_queue{vecmem_queue.queue()};
    TRACCC_INFO("Running on device: " << vecmem_queue.device_name());

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::sycl::host_memory_resource sycl_host_mr{vecmem_queue};
    vecmem::sycl::shared_memory_resource shared_mr{vecmem_queue};
    vecmem::sycl::device_memory_resource device_mr{vecmem_queue};
    traccc::memory_resource mr{device_mr, &sycl_host_mr};

    // Copy object for asynchronous data transfers.
    vecmem::sycl::async_copy copy{vecmem_queue};

    // Performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{},
        logger().clone("SeedingPerformanceWriter"));

    // Output stats
    uint64_t n_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_sycl = 0;

    /*****************************
     * Build a geometry
     *****************************/

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector host_det;
    traccc::io::read_detector(host_det, host_mr, detector_opts.detector_file,
                              detector_opts.material_file,
                              detector_opts.grid_file);

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(host_det, device_mr, copy);
    vecmem_queue.synchronize();

    const traccc::vector3 field_vec(seeding_opts);
    const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
    const auto device_field =
        traccc::sycl::make_magnetic_field(host_field, traccc_queue);

    // Seeding algorithm
    const traccc::seedfinder_config seedfinder_config(seeding_opts);
    const traccc::seedfilter_config seedfilter_config(seeding_opts);
    const traccc::spacepoint_grid_config spacepoint_grid_config(seeding_opts);
    traccc::host::seeding_algorithm sa(
        seedfinder_config, spacepoint_grid_config, seedfilter_config, host_mr,
        logger().clone("HostSeedingAlg"));
    const traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(
        track_params_estimation_config, host_mr,
        logger().clone("HostTrackParEstAlg"));

    traccc::sycl::triplet_seeding_algorithm sa_sycl{
        seedfinder_config,
        spacepoint_grid_config,
        seedfilter_config,
        mr,
        copy,
        traccc_queue,
        logger().clone("SyclSeedingAlg")};
    traccc::sycl::seed_parameter_estimation_algorithm tp_sycl{
        track_params_estimation_config, mr, copy, traccc_queue,
        logger().clone("SyclTrackParEstAlg")};

    traccc::performance::timing_info elapsedTimes;

    // Loop over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Instantiate host containers/collections
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::edm::spacepoint_collection::host spacepoints_per_event{host_mr};
        traccc::host::seeding_algorithm::output_type seeds{host_mr};
        traccc::host::track_params_estimation::output_type params{&host_mr};

        // Instantiate sycl containers/collections
        traccc::edm::seed_collection::buffer seeds_sycl_buffer;
        traccc::bound_track_parameters_collection_types::buffer
            params_sycl_buffer(0, *mr.host);

        {  // Start measuring wall time
            traccc::performance::timer wall_t("Wall time", elapsedTimes);

            /*-----------------
              hit file reading
              -----------------*/

            {
                traccc::performance::timer t("Hit reading  (cpu)",
                                             elapsedTimes);
                // Read the hits from the relevant event file
                traccc::io::read_spacepoints(
                    spacepoints_per_event, measurements_per_event, event,
                    input_opts.directory,
                    (input_opts.use_acts_geom_source ? &host_det : nullptr),
                    nullptr, nullptr, input_opts.format);

            }  // stop measuring hit reading timer

            /*----------------------------
                 Seeding algorithm
              ----------------------------*/

            /// SYCL

            // Copy the measurements and spacepoint and module data to the
            // device.
            traccc::edm::measurement_collection::buffer
                measurements_sycl_buffer(
                    static_cast<unsigned int>(measurements_per_event.size()),
                    mr.main);
            copy(vecmem::get_data(measurements_per_event),
                 measurements_sycl_buffer)
                ->wait();
            traccc::edm::spacepoint_collection::buffer spacepoints_sycl_buffer(
                static_cast<unsigned int>(spacepoints_per_event.size()),
                mr.main);
            copy(vecmem::get_data(spacepoints_per_event),
                 spacepoints_sycl_buffer)
                ->wait();

            {
                traccc::performance::timer t("Seeding (sycl)", elapsedTimes);
                // Reconstruct the spacepoints into seeds.
                seeds_sycl_buffer = sa_sycl(spacepoints_sycl_buffer);
            }  // stop measuring seeding sycl timer

            // CPU

            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Seeding  (cpu)", elapsedTimes);
                seeds = sa(vecmem::get_data(spacepoints_per_event));
            }  // stop measuring seeding cpu timer

            /*----------------------------
              Track params estimation
              ----------------------------*/

            // SYCL

            {
                traccc::performance::timer t("Track params (sycl)",
                                             elapsedTimes);
                params_sycl_buffer =
                    tp_sycl(device_field, measurements_sycl_buffer,
                            spacepoints_sycl_buffer, seeds_sycl_buffer);
            }  // stop measuring track params sycl timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track params  (cpu)",
                                             elapsedTimes);
                params = tp(vecmem::get_data(measurements_per_event),
                            vecmem::get_data(spacepoints_per_event),
                            vecmem::get_data(seeds), field_vec);
            }  // stop measuring track params cpu timer

        }  // Stop measuring wall time

        /*----------------------------------
          compare seeds from cpu and sycl
          ----------------------------------*/

        // Copy the seeds to the host for comparison.
        traccc::edm::seed_collection::host seeds_sycl{host_mr};
        traccc::bound_track_parameters_collection_types::host params_sycl{
            &host_mr};
        copy(seeds_sycl_buffer, seeds_sycl)->wait();
        copy(params_sycl_buffer, params_sycl)->wait();

        if (accelerator_opts.compare_with_cpu) {
            // Show which event we are currently presenting the results for.
            TRACCC_INFO("===>>> Event " << event << " <<<===");

            // Compare the seeds made on the host and on the device
            traccc::soa_comparator<traccc::edm::seed_collection> compare_seeds{
                "seeds", traccc::details::comparator_factory<
                             traccc::edm::seed_collection::const_device::
                                 const_proxy_type>{
                             vecmem::get_data(spacepoints_per_event),
                             vecmem::get_data(spacepoints_per_event)}};
            compare_seeds(vecmem::get_data(seeds),
                          vecmem::get_data(seeds_sycl));

            // Compare the track parameters made on the host and on the device.
            traccc::collection_comparator<traccc::bound_track_parameters<>>
                compare_track_parameters{"track parameters"};
            compare_track_parameters(vecmem::get_data(params),
                                     vecmem::get_data(params_sycl));
        }

        /*----------------
             Statistics
          ---------------*/

        n_spacepoints += spacepoints_per_event.size();
        n_seeds_sycl += seeds_sycl.size();
        n_seeds += seeds.size();

        /*------------
          Writer
          ------------*/

        if (performance_opts.run) {

            traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                        input_opts.use_acts_geom_source,
                                        &host_det, input_opts.format, false);

            sd_performance_writer.write(
                vecmem::get_data(seeds_sycl),
                vecmem::get_data(spacepoints_per_event),
                vecmem::get_data(measurements_per_event), evt_data);
        }
    }

    if (performance_opts.run) {
        sd_performance_writer.finalize();
    }

    TRACCC_INFO("==> Statistics ... ");
    TRACCC_INFO("- read    " << n_spacepoints << " spacepoints");
    TRACCC_INFO("- created  (cpu)  " << n_seeds << " seeds");
    TRACCC_INFO("- created (sycl) " << n_seeds_sycl << " seeds");
    TRACCC_INFO("==>Elapsed times... " << elapsedTimes);

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSeedingSycl", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_seeding seeding_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::program_options program_opts{
        "Full Tracking Chain Using SYCL (without clusterization)",
        {detector_opts, bfield_opts, input_opts, seeding_opts, performance_opts,
         accelerator_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(detector_opts, bfield_opts, seeding_opts, input_opts,
                   performance_opts, accelerator_opts, logger->clone());
}
