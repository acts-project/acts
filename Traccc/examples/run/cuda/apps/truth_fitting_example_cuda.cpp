/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/accelerator.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_fitting.hpp"
#include "traccc/options/track_propagation.hpp"
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
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>

using namespace traccc;

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> ilogger = traccc::getDefaultLogger(
        "TracccExampleTruthFittingCuda", traccc::Logging::Level::INFO);
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::program_options program_opts{
        "Truth Track Fitting Using CUDA",
        {detector_opts, bfield_opts, input_opts, propagation_opts,
         performance_opts, accelerator_opts},
        argc,
        argv,
        logger().cloneWithSuffix("Options")};

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::host_memory_resource cuda_host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource device_mr;
    traccc::memory_resource mr{device_mr, &cuda_host_mr};

    // Performance writer
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{},
        logger().clone("FittingPerformanceWriter"));

    // Output Stats
    std::size_t n_fitted_tracks = 0;
    std::size_t n_fitted_tracks_cuda = 0;

    /*****************************
     * Build a geometry
     *****************************/

    // B field value
    const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
    const auto device_field = traccc::cuda::make_magnetic_field(host_field);

    // Read the detector
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

    /// Standard deviations for seed track parameters
    static constexpr std::array<scalar, e_bound_size> stddevs = {
        0.03f * traccc::unit<scalar>::mm,
        0.03f * traccc::unit<scalar>::mm,
        0.017f,
        0.017f,
        0.001f / traccc::unit<scalar>::GeV,
        1.f * traccc::unit<scalar>::ns};

    // Fitting algorithm object
    traccc::fitting_config fit_cfg(fitting_opts);
    fit_cfg.propagation = propagation_opts;

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
                    truth_track_candidates, truth_measurements, sg, host_mr);
            });
        truth_track_candidates.measurements =
            vecmem::get_data(truth_measurements);

        // track candidates buffer
        auto truth_measurements_buffer =
            async_copy.to(vecmem::get_data(truth_measurements), mr.main,
                          mr.host, vecmem::copy::type::host_to_device);
        traccc::edm::track_container<traccc::default_algebra>::buffer
            truth_track_candidates_buffer{
                async_copy.to(vecmem::get_data(truth_track_candidates.tracks),
                              mr.main, mr.host,
                              vecmem::copy::type::host_to_device),
                {},
                truth_measurements_buffer};

        // Instantiate cuda containers/collections
        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_states_cuda_buffer;

        {
            traccc::performance::timer t("Track fitting  (cuda)", elapsedTimes);

            // Run fitting
            track_states_cuda_buffer = device_fitting(
                detector_buffer, device_field, truth_track_candidates_buffer);
        }

        traccc::edm::track_container<traccc::default_algebra>::host
            track_states_cuda{host_mr};
        async_copy(track_states_cuda_buffer.tracks, track_states_cuda.tracks,
                   vecmem::copy::type::device_to_host)
            ->wait();
        async_copy(track_states_cuda_buffer.states, track_states_cuda.states,
                   vecmem::copy::type::device_to_host)
            ->wait();

        // CPU container(s)
        traccc::host::kalman_fitting_algorithm::output_type track_states{
            host_mr};

        if (accelerator_opts.compare_with_cpu) {

            {
                traccc::performance::timer t("Track fitting  (cpu)",
                                             elapsedTimes);

                // Run fitting
                track_states = host_fitting(
                    polymorphic_detector, host_field,
                    traccc::edm::track_container<traccc::default_algebra>::
                        const_data(truth_track_candidates));
            }
        }

        if (accelerator_opts.compare_with_cpu) {
            // Show which event we are currently presenting the results for.
            std::cout << "===>>> Event " << event << " <<<===" << std::endl;

            // Compare the track parameters made on the host and on the device.
            traccc::soa_comparator<
                traccc::edm::track_collection<traccc::default_algebra>>
                compare_track_fits{
                    "track fits",
                    traccc::details::comparator_factory<
                        traccc::edm::track_collection<traccc::default_algebra>::
                            const_device::const_proxy_type>{
                        truth_track_candidates.measurements,
                        truth_track_candidates.measurements,
                        vecmem::get_data(track_states.states),
                        vecmem::get_data(track_states_cuda.states)}};
            compare_track_fits(vecmem::get_data(track_states.tracks),
                               vecmem::get_data(track_states_cuda.tracks));
        }

        // Statistics
        n_fitted_tracks += track_states.tracks.size();
        n_fitted_tracks_cuda += track_states_cuda.tracks.size();

        if (performance_opts.run) {
            for (unsigned int i = 0; i < track_states_cuda.tracks.size(); i++) {
                host_detector_visitor<detector_type_list>(
                    polymorphic_detector,
                    [&]<typename detector_traits_t>(
                        const typename detector_traits_t::host& det) {
                        fit_performance_writer.write(
                            track_states_cuda.tracks.at(i),
                            track_states_cuda.states, truth_measurements, det,
                            evt_data);
                    });
            }
        }
    }

    if (performance_opts.run) {
        fit_performance_writer.finalize();
    }

    std::cout << "==> Statistics ... " << std::endl;
    std::cout << "- created (cuda) " << n_fitted_tracks_cuda << " fitted tracks"
              << std::endl;
    std::cout << "- created  (cpu) " << n_fitted_tracks << " fitted tracks"
              << std::endl;
    std::cout << "==>Elapsed times...\n" << elapsedTimes << std::endl;

    return EXIT_SUCCESS;
}
