/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/efficiency/nseed_performance_writer.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/efficiency/track_filter.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/examples/print_fitted_tracks_statistics.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/accelerator.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
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

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/options/track_gbts_seeding.hpp"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

using namespace traccc;

int seq_run(const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::track_gbts_seeding& seeding_gbts_opts,
            const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            [[maybe_unused]] const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            const traccc::opts::truth_finding& truth_finding_opts,
            const traccc::opts::seed_matching& seed_matching_opts,
            const traccc::opts::track_matching& track_matching_opts,
            std::unique_ptr<const traccc::Logger> ilogger, bool usingGBTS) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::host_memory_resource cuda_host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr{device_mr, &cuda_host_mr};

    // Performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .seed_truth_config = seed_matching_opts},
        logger().clone("SeedingPerformanceWriter"));
    traccc::finding_performance_writer find_performance_writer(
        traccc::finding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .track_truth_config = track_matching_opts,
            .require_fit = true},
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
    uint64_t n_seeds_cuda = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_found_tracks_cuda = 0;
    uint64_t n_fitted_tracks = 0;
    uint64_t n_fitted_tracks_cuda = 0;

    /*****************************
     * Build a geometry
     *****************************/

    traccc::detector_design_description::host host_det_descr{host_mr};
    traccc::detector_conditions_description::host host_det_cond{host_mr};
    traccc::io::read_detector_description(
        host_det_descr, host_det_cond, detector_opts.detector_file,
        detector_opts.digitization_file, detector_opts.conditions_file,
        traccc::data_format::json);

    // B field value
    const traccc::vector3 field_vec(seeding_opts);
    const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
    const auto device_field = traccc::cuda::make_magnetic_field(
        host_field,
        (accelerator_opts.use_gpu_texture_memory
             ? traccc::cuda::magnetic_field_storage::texture_memory
             : traccc::cuda::magnetic_field_storage::global_memory));

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector host_det;
    traccc::io::read_detector(host_det, mng_mr, detector_opts.detector_file,
                              detector_opts.material_file,
                              detector_opts.grid_file);

    // Copy objects
    vecmem::copy host_copy;
    vecmem::cuda::copy copy;

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(host_det, mng_mr, host_copy);

    // GBTS seeding configuration
    const traccc::gbts_seedfinder_config gbts_config(seeding_gbts_opts);
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

    traccc::cuda::stream stream;

    vecmem::cuda::async_copy async_copy{stream.cudaStream()};

    traccc::cuda::triplet_seeding_algorithm sa_cuda{
        seedfinder_config,
        spacepoint_grid_config,
        seedfilter_config,
        mr,
        async_copy,
        stream,
        logger().clone("CudaSeedingAlg")};
    traccc::cuda::gbts_seeding_algorithm gbts_sa_cuda(
        gbts_config, mr, copy, stream, logger().clone("CudaGbtsSeedingAlg"));
    traccc::cuda::seed_parameter_estimation_algorithm tp_cuda{
        track_params_estimation_config, mr, async_copy, stream,
        logger().clone("CudaTrackParEstAlg")};

    // Propagation configuration
    detray::propagation::config propagation_config(propagation_opts);

    // Finding algorithm configuration
    traccc::finding_config cfg(finding_opts);
    cfg.propagation = propagation_config;

    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(
        cfg, host_mr, logger().clone("HostFindingAlg"));
    traccc::cuda::combinatorial_kalman_filter_algorithm device_finding(
        cfg, mr, async_copy, stream, logger().clone("CudaFindingAlg"));

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

        traccc::edm::seed_collection::buffer seeds_cuda_buffer;
        traccc::bound_track_parameters_collection_types::buffer
            params_cuda_buffer(0, *mr.host);

        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_candidates_cuda_buffer;

        {  // Start measuring wall time
            traccc::performance::timer wall_t("Wall time", elapsedTimes);

            /*-----------------
            hit file reading
            -----------------*/
            {
                traccc::performance::timer t("Hit reading  (cpu)",
                                             elapsedTimes);
                // Read the hits and measurements from the relevant event files
                traccc::io::read_spacepoints(
                    spacepoints_per_event, measurements_per_event, event,
                    input_opts.directory,
                    (input_opts.use_acts_geom_source ? &host_det : nullptr),
                    &host_det_descr, &host_det_cond, input_opts.format);

            }  // stop measuring hit reading timer

            /*----------------------------
                Seeding algorithm
            ----------------------------*/

            /// CUDA

            // Copy the spacepoint and module data to the device.
            traccc::edm::spacepoint_collection::buffer spacepoints_cuda_buffer(
                static_cast<unsigned int>(spacepoints_per_event.size()),
                mr.main);
            async_copy.setup(spacepoints_cuda_buffer)->wait();
            async_copy(vecmem::get_data(spacepoints_per_event),
                       spacepoints_cuda_buffer)
                ->wait();

            traccc::edm::measurement_collection::buffer
                measurements_cuda_buffer(
                    static_cast<unsigned int>(measurements_per_event.size()),
                    mr.main);
            async_copy.setup(measurements_cuda_buffer)->wait();
            async_copy(vecmem::get_data(measurements_per_event),
                       measurements_cuda_buffer)
                ->wait();

            {
                traccc::performance::timer t("Seeding (cuda)", elapsedTimes);
                // Reconstruct the spacepoints into seeds.
                if (usingGBTS) {
                    seeds_cuda_buffer = gbts_sa_cuda(spacepoints_cuda_buffer,
                                                     measurements_cuda_buffer);
                } else {
                    seeds_cuda_buffer = sa_cuda(spacepoints_cuda_buffer);
                }
                stream.synchronize();
            }  // stop measuring seeding cuda timer

            // CPU

            if (accelerator_opts.compare_with_cpu) {
                {
                    traccc::performance::timer t("Seeding  (cpu)",
                                                 elapsedTimes);
                    seeds = sa(vecmem::get_data(spacepoints_per_event));
                }
            }  // stop measuring seeding cpu timer

            /*----------------------------
               Track params estimation
            ----------------------------*/

            // CUDA
            {
                traccc::performance::timer t("Track params (cuda)",
                                             elapsedTimes);
                params_cuda_buffer =
                    tp_cuda(device_field, measurements_cuda_buffer,
                            spacepoints_cuda_buffer, seeds_cuda_buffer);
                stream.synchronize();
            }  // stop measuring track params cuda timer

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
                traccc::performance::timer t("Track finding with CKF (cuda)",
                                             elapsedTimes);
                track_candidates_cuda_buffer = device_finding(
                    detector_buffer, device_field, measurements_cuda_buffer,
                    params_cuda_buffer);
            }

            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track finding with CKF (cpu)",
                                             elapsedTimes);
                track_candidates =
                    host_finding(host_det, host_field,
                                 vecmem::get_data(measurements_per_event),
                                 vecmem::get_data(params));
            }
        }  // Stop measuring wall time

        /*----------------------------------
          compare seeds from cpu and cuda
          ----------------------------------*/

        // Copy the seeds to the host for comparisons
        traccc::edm::seed_collection::host seeds_cuda{host_mr};
        traccc::bound_track_parameters_collection_types::host params_cuda;
        async_copy(seeds_cuda_buffer, seeds_cuda)->wait();
        async_copy(params_cuda_buffer, params_cuda)->wait();

        // Copy track candidates from device to host
        traccc::edm::track_container<traccc::default_algebra>::host
            track_candidates_cuda{host_mr,
                                  vecmem::get_data(measurements_per_event)};
        async_copy(track_candidates_cuda_buffer.tracks,
                   track_candidates_cuda.tracks)
            ->wait();
        async_copy(track_candidates_cuda_buffer.states,
                   track_candidates_cuda.states)
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
                          vecmem::get_data(seeds_cuda));

            // Compare the track parameters made on the host and on the device.
            traccc::collection_comparator<traccc::bound_track_parameters<>>
                compare_track_parameters{"track parameters"};
            compare_track_parameters(vecmem::get_data(params),
                                     vecmem::get_data(params_cuda));
        }

        /*----------------
             Statistics
          ---------------*/

        details::print_fitted_tracks_statistics(track_candidates_cuda,
                                                logger());
        n_spacepoints += spacepoints_per_event.size();
        n_seeds_cuda += seeds_cuda.size();
        n_seeds += seeds.size();
        n_found_tracks_cuda += track_candidates_cuda.tracks.size();
        n_found_tracks += track_candidates.tracks.size();
        n_fitted_tracks_cuda += track_candidates_cuda.tracks.size();
        n_fitted_tracks += track_candidates.tracks.size();

        /*------------
          Writer
          ------------*/

        if (performance_opts.run) {

            traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                        input_opts.use_acts_geom_source,
                                        &host_det, input_opts.format, false);

            sd_performance_writer.write(
                vecmem::get_data(seeds_cuda),
                vecmem::get_data(spacepoints_per_event),
                vecmem::get_data(measurements_per_event), evt_data);

            std::cout << track_candidates_cuda.tracks.size() << ", "
                      << track_candidates_cuda.states.size() << std::endl;

            find_performance_writer.write(
                traccc::edm::track_container<
                    traccc::default_algebra>::const_data(track_candidates_cuda),
                evt_data);

            for (unsigned int i = 0; i < track_candidates_cuda.tracks.size();
                 i++) {
                host_detector_visitor<detector_type_list>(
                    host_det, [&]<typename detector_traits_t>(
                                  const typename detector_traits_t::host& det) {
                        fit_performance_writer.write(
                            track_candidates_cuda.tracks.at(i),
                            track_candidates_cuda.states,
                            measurements_per_event, det, evt_data);
                    });
            }
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
    std::cout << "- created (cuda)  " << n_seeds_cuda << " seeds" << std::endl;
    std::cout << "- created  (cpu) " << n_found_tracks << " found tracks"
              << std::endl;
    std::cout << "- created (cuda) " << n_found_tracks_cuda << " found tracks"
              << std::endl;
    std::cout << "- created  (cpu) " << n_fitted_tracks << " fitted tracks"
              << std::endl;
    std::cout << "- created (cuda) " << n_fitted_tracks_cuda << " fitted tracks"
              << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSeedingCuda", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_seeding seeding_opts;
    traccc::opts::track_gbts_seeding seeding_gbts_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::truth_finding truth_finding_opts;
    traccc::opts::seed_matching seed_matching_opts;
    traccc::opts::track_matching track_matching_opts;
    traccc::opts::program_options program_opts{
        "Full Tracking Chain Using CUDA (without clusterization)",
        {detector_opts, bfield_opts, input_opts, seeding_opts,
         seeding_gbts_opts, finding_opts, propagation_opts, fitting_opts,
         performance_opts, accelerator_opts, truth_finding_opts,
         seed_matching_opts, track_matching_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(seeding_opts, seeding_gbts_opts, finding_opts,
                   propagation_opts, fitting_opts, input_opts, detector_opts,
                   bfield_opts, performance_opts, accelerator_opts,
                   truth_finding_opts, seed_matching_opts, track_matching_opts,
                   logger->clone(), seeding_gbts_opts.useGBTS);
}
