/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/alpaka/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/alpaka/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/alpaka/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/alpaka/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/alpaka/utils/queue.hpp"
#include "traccc/alpaka/utils/vecmem_objects.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/efficiency/nseed_performance_writer.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/efficiency/track_filter.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/accelerator.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/seed_matching.hpp"
#include "traccc/options/track_finding.hpp"
#include "traccc/options/track_fitting.hpp"
#include "traccc/options/track_matching.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/options/track_seeding.hpp"
#include "traccc/options/truth_finding.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/soa_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/utils/propagation.hpp"

// System include(s).
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>

using namespace traccc;

int seq_run(const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            const traccc::opts::truth_finding& truth_finding_opts,
            const traccc::opts::seed_matching& seed_matching_opts,
            const traccc::opts::track_matching& track_matching_opts,
            [[maybe_unused]] std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resources used by the application.
    traccc::alpaka::queue queue;
    traccc::alpaka::vecmem_objects vo(queue);

    vecmem::memory_resource& host_mr = vo.host_mr();
    vecmem::memory_resource& device_mr = vo.device_mr();
    vecmem::memory_resource& mng_mr = vo.shared_mr();
    traccc::memory_resource mr{device_mr, &host_mr};

    // Performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .seed_truth_config = seed_matching_opts},
        logger().clone("SeedingPerformanceWriter"));
    traccc::finding_performance_writer find_performance_writer(
        traccc::finding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .track_truth_config = track_matching_opts},
        logger().clone("FindingPerformanceWriter"));
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{},
        logger().clone("FittingPerformanceWriter"));

    traccc::nseed_performance_writer nsd_performance_writer(
        "nseed_performance_",
        std::make_unique<traccc::simple_charged_eta_pt_cut>(
            2.7f, 1.f * traccc::unit<traccc::scalar>::GeV),
        std::make_unique<traccc::stepped_percentage>(0.6f));

    if (performance_opts.run) {
        nsd_performance_writer.initialize();
    }

    // Output stats
    uint64_t n_spacepoints = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_alpaka = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_found_tracks_alpaka = 0;
    uint64_t n_fitted_tracks = 0;
    uint64_t n_fitted_tracks_alpaka = 0;

    /*****************************
     * Build a geometry
     *****************************/

    // B field value and its type
    const auto field = traccc::details::make_magnetic_field(bfield_opts);
    const traccc::vector3 field_vec(seeding_opts);

    // Detector view object
    traccc::host_detector host_det;
    traccc::io::read_detector(host_det, mng_mr, detector_opts.detector_file,
                              detector_opts.material_file,
                              detector_opts.grid_file);

    // Copy objects
    vecmem::copy host_copy;
    vecmem::copy& copy = vo.copy();
    vecmem::copy& async_copy = vo.async_copy();

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(host_det, mng_mr, copy);

    // Seeding algorithms
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

    // Alpaka Algorithms
    traccc::alpaka::triplet_seeding_algorithm sa_alpaka{
        seedfinder_config,
        spacepoint_grid_config,
        seedfilter_config,
        mr,
        async_copy,
        queue,
        logger().clone("AlpakaSeedingAlg")};
    traccc::alpaka::seed_parameter_estimation_algorithm tp_alpaka{
        track_params_estimation_config, mr, async_copy, queue,
        logger().clone("AlpakaTrackParEstAlg")};

    // Propagation configuration
    detray::propagation::config propagation_config(propagation_opts);

    // Finding algorithm configuration
    traccc::finding_config cfg(finding_opts);
    cfg.propagation = propagation_config;

    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(
        cfg, host_mr, logger().clone("HostFindingAlg"));
    traccc::alpaka::combinatorial_kalman_filter_algorithm device_finding(
        cfg, mr, copy, queue, logger().clone("AlpakaFindingAlg"));

    // Fitting algorithm object
    traccc::fitting_config fit_cfg(fitting_opts);
    fit_cfg.propagation = propagation_config;

    traccc::host::kalman_fitting_algorithm host_fitting(
        fit_cfg, host_mr, host_copy, logger().clone("HostFittingAlg"));
    traccc::alpaka::kalman_fitting_algorithm device_fitting(
        fit_cfg, mr, copy, queue, logger().clone("AlpakaFittingAlg"));

    traccc::performance::timing_info elapsedTimes;

    // Loop over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Instantiate host containers/collections
        traccc::edm::spacepoint_collection::host spacepoints_per_event{host_mr};
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::host::seeding_algorithm::output_type seeds{host_mr};
        traccc::host::track_params_estimation::output_type params;
        traccc::edm::track_container<traccc::default_algebra>::host
            track_candidates{host_mr};
        traccc::edm::track_container<traccc::default_algebra>::host
            track_states{host_mr};

        traccc::edm::seed_collection::buffer seeds_alpaka_buffer;
        traccc::bound_track_parameters_collection_types::buffer
            params_alpaka_buffer(0, *mr.host);

        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_candidates_alpaka_buffer;

        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_states_alpaka_buffer;

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

            // Alpaka

            // Copy the spacepoint data to the device.
            traccc::edm::spacepoint_collection::buffer
                spacepoints_alpaka_buffer(
                    static_cast<unsigned int>(spacepoints_per_event.size()),
                    mr.main);
            async_copy.setup(spacepoints_alpaka_buffer)->wait();
            async_copy(vecmem::get_data(spacepoints_per_event),
                       spacepoints_alpaka_buffer)
                ->wait();

            traccc::edm::measurement_collection::buffer
                measurements_alpaka_buffer(
                    static_cast<unsigned int>(measurements_per_event.size()),
                    mr.main);
            async_copy.setup(measurements_alpaka_buffer)->wait();
            async_copy(vecmem::get_data(measurements_per_event),
                       measurements_alpaka_buffer)
                ->wait();

            {
                traccc::performance::timer t("Seeding (alpaka)", elapsedTimes);
                // Reconstruct the spacepoints into seeds.
                seeds_alpaka_buffer =
                    sa_alpaka(vecmem::get_data(spacepoints_alpaka_buffer));
                queue.synchronize();
            }

            // CPU

            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Seeding  (cpu)", elapsedTimes);
                seeds = sa(vecmem::get_data(spacepoints_per_event));
            }  // stop measuring seeding cpu timer

            /*----------------------------
            Track params estimation
            ----------------------------*/

            // Alpaka

            {
                traccc::performance::timer t("Track params (alpaka)",
                                             elapsedTimes);
                params_alpaka_buffer =
                    tp_alpaka(field, measurements_alpaka_buffer,
                              spacepoints_alpaka_buffer, seeds_alpaka_buffer);
                queue.synchronize();
            }  // stop measuring track params alpaka timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track params  (cpu)",
                                             elapsedTimes);
                params = tp(vecmem::get_data(measurements_per_event),
                            vecmem::get_data(spacepoints_per_event),
                            vecmem::get_data(seeds), field_vec);
            }  // stop measuring track params cpu timer

            /*------------------------
               Track Finding with CKF
              ------------------------*/

            {
                traccc::performance::timer t("Track finding with CKF (alpaka)",
                                             elapsedTimes);
                track_candidates_alpaka_buffer = device_finding(
                    detector_buffer, field, measurements_alpaka_buffer,
                    params_alpaka_buffer);
            }

            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track finding with CKF (cpu)",
                                             elapsedTimes);
                track_candidates = host_finding(
                    host_det, field, vecmem::get_data(measurements_per_event),
                    vecmem::get_data(params));
            }

            /*------------------------
               Track Fitting with KF
              ------------------------*/

            {
                traccc::performance::timer t("Track fitting with KF (alpaka)",
                                             elapsedTimes);

                track_states_alpaka_buffer = device_fitting(
                    detector_buffer, field, track_candidates_alpaka_buffer);
            }

            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track fitting with KF (cpu)",
                                             elapsedTimes);
                track_states = host_fitting(
                    host_det, field,
                    traccc::edm::track_container<
                        traccc::default_algebra>::const_data(track_candidates));
            }

        }  // Stop measuring wall time

        /*----------------------------------
          compare seeds from cpu and alpaka
          ----------------------------------*/

        // Copy the seeds to the host for comparisons
        traccc::edm::seed_collection::host seeds_alpaka{host_mr};
        traccc::bound_track_parameters_collection_types::host params_alpaka{
            &host_mr};
        async_copy(seeds_alpaka_buffer, seeds_alpaka)->wait();
        async_copy(params_alpaka_buffer, params_alpaka)->wait();

        // Copy track candidates from device to host
        traccc::edm::track_collection<traccc::default_algebra>::host
            track_candidates_alpaka{host_mr};
        copy(track_candidates_alpaka_buffer.tracks, track_candidates_alpaka)
            ->wait();

        // Copy track states from device to host
        traccc::edm::track_container<traccc::default_algebra>::host
            track_states_alpaka{host_mr};
        async_copy(track_states_alpaka_buffer.tracks,
                   track_states_alpaka.tracks)
            ->wait();
        async_copy(track_states_alpaka_buffer.states,
                   track_states_alpaka.states)
            ->wait();

        if (accelerator_opts.compare_with_cpu) {
            // Show which event we are currently presenting the results for.
            std::cout << "===>>> Event " << event << " <<<===" << std::endl;

            // Compare the seeds made on the host and on the device
            traccc::soa_comparator<traccc::edm::seed_collection> compare_seeds{
                "seeds", traccc::details::comparator_factory<
                             traccc::edm::seed_collection::const_device::
                                 const_proxy_type>{
                             vecmem::get_data(spacepoints_per_event),
                             vecmem::get_data(spacepoints_per_event)}};
            compare_seeds(vecmem::get_data(seeds),
                          vecmem::get_data(seeds_alpaka));

            // Compare the track parameters made on the host and on the device.
            traccc::collection_comparator<traccc::bound_track_parameters<>>
                compare_track_parameters{"track parameters"};
            compare_track_parameters(vecmem::get_data(params),
                                     vecmem::get_data(params_alpaka));

            // Compare the track candidates made on the host and on the
            // device
            traccc::soa_comparator<
                traccc::edm::track_collection<traccc::default_algebra>>
                compare_track_candidates{
                    "track candidates",
                    traccc::details::comparator_factory<
                        traccc::edm::track_collection<traccc::default_algebra>::
                            const_device::const_proxy_type>{
                        vecmem::get_data(measurements_per_event),
                        vecmem::get_data(measurements_per_event),
                        {},
                        {}}};
            compare_track_candidates(vecmem::get_data(track_candidates.tracks),
                                     vecmem::get_data(track_candidates_alpaka));
        }

        /*----------------
             Statistics
          ---------------*/

        n_spacepoints += spacepoints_per_event.size();
        n_seeds_alpaka += seeds_alpaka.size();
        n_seeds += seeds.size();
        n_found_tracks_alpaka += track_candidates_alpaka.size();
        n_found_tracks += track_candidates.tracks.size();
        n_fitted_tracks_alpaka += track_states_alpaka.tracks.size();
        n_fitted_tracks += track_states.tracks.size();

        /*------------
          Writer
          ------------*/

        if (performance_opts.run) {

            traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                        input_opts.use_acts_geom_source,
                                        &host_det, input_opts.format, false);

            sd_performance_writer.write(
                vecmem::get_data(seeds),
                vecmem::get_data(spacepoints_per_event),
                vecmem::get_data(measurements_per_event), evt_data);
        }
    }

    if (performance_opts.run) {
        sd_performance_writer.finalize();
        nsd_performance_writer.finalize();
        find_performance_writer.finalize();
        fit_performance_writer.finalize();
        std::cout << nsd_performance_writer.generate_report_str();
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- read    " << n_spacepoints << " spacepoints" << std::endl;
    std::cout << "- created  (cpu)  " << n_seeds << " seeds" << std::endl;
    std::cout << "- created (alpaka)  " << n_seeds_alpaka << " seeds"
              << std::endl;
    std::cout << "- created  (cpu) " << n_found_tracks << " found tracks"
              << std::endl;
    std::cout << "- created (alpaka) " << n_found_tracks_alpaka
              << " found tracks" << std::endl;
    std::cout << "- created  (cpu) " << n_fitted_tracks << " fitted tracks"
              << std::endl;
    std::cout << "- created (alpaka) " << n_fitted_tracks_alpaka
              << " fitted tracks" << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSeedingAlpaka", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_seeding seeding_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::truth_finding truth_finding_opts;
    traccc::opts::seed_matching seed_matching_opts;
    traccc::opts::track_matching track_matching_opts;
    traccc::opts::program_options program_opts{
        "Full Tracking Chain Using Alpaka (without clusterization)",
        {detector_opts, bfield_opts, input_opts, seeding_opts, finding_opts,
         propagation_opts, fitting_opts, performance_opts, accelerator_opts,
         truth_finding_opts, seed_matching_opts, track_matching_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(seeding_opts, finding_opts, propagation_opts, fitting_opts,
                   input_opts, detector_opts, bfield_opts, performance_opts,
                   accelerator_opts, truth_finding_opts, seed_matching_opts,
                   track_matching_opts, logger->clone());
}
