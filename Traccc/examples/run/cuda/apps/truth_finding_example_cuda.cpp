/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/accelerator.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_finding.hpp"
#include "traccc/options/track_fitting.hpp"
#include "traccc/options/track_matching.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/options/truth_finding.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/soa_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/utils/propagation.hpp"
#include "traccc/utils/seed_generator.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>

using namespace traccc;

int seq_run(const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            const traccc::opts::truth_finding& truth_finding_opts,
            const traccc::opts::track_matching& track_matching_opts,
            std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::host_memory_resource cuda_host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr{device_mr, &cuda_host_mr};

    // Performance writer
    traccc::finding_performance_writer find_performance_writer(
        traccc::finding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .track_truth_config = track_matching_opts},
        logger().clone("FindingPerformanceWriter"));
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{},
        logger().clone("FittingPerformanceWriter"));

    // Output Stats
    uint64_t n_found_tracks = 0;
    uint64_t n_found_tracks_cuda = 0;
    uint64_t n_fitted_tracks = 0;
    uint64_t n_fitted_tracks_cuda = 0;

    /*****************************
     * Build a geometry
     *****************************/

    // B field value
    const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
    const auto device_field = traccc::cuda::make_magnetic_field(host_field);

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector polymorphic_detector;
    traccc::io::read_detector(
        polymorphic_detector, mng_mr, detector_opts.detector_file,
        detector_opts.material_file, detector_opts.grid_file);

    /*****************************
     * Do the reconstruction
     *****************************/

    // Stream object
    traccc::cuda::stream stream;

    // Copy object
    vecmem::copy host_copy;
    vecmem::cuda::async_copy async_copy{stream.cudaStream()};

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(polymorphic_detector, device_mr,
                                          async_copy);

    // Standard deviations for seed track parameters
    static constexpr std::array<traccc::scalar, traccc::e_bound_size> stddevs =
        {1e-4f * traccc::unit<traccc::scalar>::mm,
         1e-4f * traccc::unit<traccc::scalar>::mm,
         1e-3f,
         1e-3f,
         1e-4f / traccc::unit<traccc::scalar>::GeV,
         1e-4f * traccc::unit<traccc::scalar>::ns};

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

    // Fitting algorithm object
    traccc::fitting_config fit_cfg(fitting_opts);
    fit_cfg.propagation = propagation_config;

    traccc::host::kalman_fitting_algorithm host_fitting(
        fit_cfg, host_mr, host_copy, logger().clone("HostFittingAlg"));
    traccc::cuda::kalman_fitting_algorithm device_fitting(
        fit_cfg, mr, async_copy, stream, logger().clone("CudaFittingAlg"));

    traccc::performance::timing_info elapsedTimes;

    // Iterate over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Truth Track Candidates
        traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                    input_opts.use_acts_geom_source,
                                    &polymorphic_detector, input_opts.format,
                                    false);

        traccc::edm::measurement_collection::host truth_measurements{host_mr};
        traccc::edm::track_container<traccc::default_algebra>::host
            truth_track_candidates{host_mr};

        host_detector_visitor<detector_type_list>(
            polymorphic_detector,
            [&]<typename detector_traits_t>(
                const typename detector_traits_t::host& det) {
                // Seed generator
                traccc::seed_generator<typename detector_traits_t::host> sg(
                    det, stddevs);
                evt_data.generate_truth_candidates(
                    truth_track_candidates, truth_measurements, sg, host_mr,
                    truth_finding_opts.m_pT_min);
            });
        truth_track_candidates.measurements =
            vecmem::get_data(truth_measurements);

        // Prepare truth seeds
        traccc::bound_track_parameters_collection_types::host seeds(mr.host);
        const std::size_t n_tracks = truth_track_candidates.tracks.size();
        for (std::size_t i_trk = 0; i_trk < n_tracks; i_trk++) {
            seeds.push_back(truth_track_candidates.tracks.at(i_trk).params());
        }

        std::cout << "Number of seeds: " << seeds.size() << std::endl;

        traccc::bound_track_parameters_collection_types::buffer seeds_buffer{
            static_cast<unsigned int>(seeds.size()), mr.main};
        async_copy.setup(seeds_buffer)->wait();
        async_copy(vecmem::get_data(seeds), seeds_buffer,
                   vecmem::copy::type::host_to_device)
            ->wait();

        // Read measurements
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::io::read_measurements(
            measurements_per_event, event, input_opts.directory,
            (input_opts.use_acts_geom_source ? &polymorphic_detector : nullptr),
            nullptr, nullptr, input_opts.format);

        traccc::edm::measurement_collection::buffer measurements_cuda_buffer(
            static_cast<unsigned int>(measurements_per_event.size()), mr.main);
        async_copy.setup(measurements_cuda_buffer)->wait();
        async_copy(vecmem::get_data(measurements_per_event),
                   measurements_cuda_buffer)
            ->wait();

        // Instantiate output cuda containers/collections
        traccc::cuda::combinatorial_kalman_filter_algorithm::output_type
            track_candidates_cuda_buffer;

        {
            traccc::performance::timer t("Track finding  (cuda)", elapsedTimes);

            // Run finding
            track_candidates_cuda_buffer =
                device_finding(detector_buffer, device_field,
                               measurements_cuda_buffer, seeds_buffer);
        }

        traccc::edm::track_container<traccc::default_algebra>::host
            track_candidates_cuda{host_mr,
                                  vecmem::get_data(measurements_per_event)};
        async_copy(track_candidates_cuda_buffer.tracks,
                   track_candidates_cuda.tracks,
                   vecmem::copy::type::device_to_host)
            ->wait();

        // Instantiate cuda containers/collections
        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_states_cuda_buffer;

        {
            traccc::performance::timer t("Track fitting  (cuda)", elapsedTimes);

            // Run fitting
            track_states_cuda_buffer = device_fitting(
                detector_buffer, device_field, track_candidates_cuda_buffer);
        }
        traccc::edm::track_container<traccc::default_algebra>::host
            track_states_cuda{host_mr};
        async_copy(track_states_cuda_buffer.tracks, track_states_cuda.tracks,
                   vecmem::copy::type::device_to_host)
            ->wait();
        async_copy(track_states_cuda_buffer.states, track_states_cuda.states,
                   vecmem::copy::type::device_to_host)
            ->wait();

        // CPU containers
        traccc::host::combinatorial_kalman_filter_algorithm::output_type
            track_candidates{host_mr};
        traccc::host::kalman_fitting_algorithm::output_type track_states{
            host_mr};

        if (accelerator_opts.compare_with_cpu) {

            {
                traccc::performance::timer t("Track finding  (cpu)",
                                             elapsedTimes);

                // Run finding
                track_candidates =
                    host_finding(polymorphic_detector, host_field,
                                 vecmem::get_data(measurements_per_event),
                                 vecmem::get_data(seeds));
            }

            {
                traccc::performance::timer t("Track fitting  (cpu)",
                                             elapsedTimes);

                // Run fitting
                track_states = host_fitting(
                    polymorphic_detector, host_field,
                    traccc::edm::track_container<
                        traccc::default_algebra>::const_data(track_candidates));
            }
        }

        if (accelerator_opts.compare_with_cpu) {

            // Show which event we are currently presenting the results for.
            TRACCC_INFO("===>>> Event " << event << " <<<===");

            // Compare the track parameters made on the host and on the device.
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
            compare_track_candidates(
                vecmem::get_data(track_candidates.tracks),
                vecmem::get_data(track_candidates_cuda.tracks));
        }

        /// Statistics
        n_found_tracks += track_candidates.tracks.size();
        n_fitted_tracks += track_states.tracks.size();
        n_found_tracks_cuda += track_candidates_cuda.tracks.size();
        n_fitted_tracks_cuda += track_states_cuda.tracks.size();

        if (performance_opts.run) {
            find_performance_writer.write(
                traccc::edm::track_container<
                    traccc::default_algebra>::const_data(track_candidates_cuda),

                evt_data);

            for (unsigned int i = 0; i < track_states_cuda.tracks.size(); i++) {
                host_detector_visitor<detector_type_list>(
                    polymorphic_detector,
                    [&]<typename detector_traits_t>(
                        const typename detector_traits_t::host& det) {
                        fit_performance_writer.write(
                            track_states_cuda.tracks.at(i),
                            track_states_cuda.states, measurements_per_event,
                            det, evt_data);
                    });
            }
        }
    }

    if (performance_opts.run) {
        find_performance_writer.finalize();
        fit_performance_writer.finalize();
    }

    TRACCC_INFO("==> Statistics ... ");
    TRACCC_INFO("- created (cuda) " << n_found_tracks_cuda << " found tracks");
    TRACCC_INFO("- created (cuda) " << n_fitted_tracks_cuda
                                    << " fitted tracks");
    TRACCC_INFO("- created  (cpu) " << n_found_tracks << " found tracks");
    TRACCC_INFO("- created  (cpu) " << n_fitted_tracks << " fitted tracks");
    TRACCC_INFO("==>Elapsed times... " << elapsedTimes);

    return 1;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleTruthFindingCuda", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::truth_finding truth_finding_config;
    traccc::opts::track_matching track_matching_opts;
    traccc::opts::program_options program_opts{
        "Truth Track Finding Using CUDA",
        {detector_opts, bfield_opts, input_opts, finding_opts, propagation_opts,
         fitting_opts, performance_opts, accelerator_opts, truth_finding_config,
         track_matching_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(finding_opts, propagation_opts, fitting_opts, input_opts,
                   detector_opts, bfield_opts, performance_opts,
                   accelerator_opts, truth_finding_config, track_matching_opts,
                   logger->clone());
}
