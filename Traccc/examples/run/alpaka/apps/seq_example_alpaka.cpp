/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/alpaka/clusterization/clusterization_algorithm.hpp"
#include "traccc/alpaka/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/alpaka/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/alpaka/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/alpaka/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/alpaka/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/alpaka/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/alpaka/utils/queue.hpp"
#include "traccc/alpaka/utils/vecmem_objects.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/accelerator.hpp"
#include "traccc/options/clusterization.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_finding.hpp"
#include "traccc/options/track_fitting.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/options/track_seeding.hpp"
#include "traccc/performance/collection_comparator.hpp"
#include "traccc/performance/container_comparator.hpp"
#include "traccc/performance/soa_comparator.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/utils/propagation.hpp"

// System include(s).
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>

int seq_run(const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::clusterization& clusterization_opts,
            const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::accelerator& accelerator_opts,
            std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resources used by the application.
    traccc::alpaka::queue queue;
    traccc::alpaka::vecmem_objects vo(queue);

    vecmem::memory_resource& host_mr = vo.host_mr();
    vecmem::memory_resource& device_mr = vo.device_mr();
    traccc::memory_resource mr{device_mr, &host_mr};

    // Host copy object
    vecmem::copy host_copy;

    // Device types used.
    vecmem::copy& copy = vo.async_copy();

    // Construct the detector description object.
    traccc::detector_design_description::host host_det_descr{host_mr};
    traccc::detector_conditions_description::host host_det_cond{host_mr};
    traccc::io::read_detector_description(
        host_det_descr, host_det_cond, detector_opts.detector_file,
        detector_opts.digitization_file, detector_opts.conditions_file,
        traccc::data_format::json);
    traccc::detector_design_description::data host_det_descr_data{
        vecmem::get_data(host_det_descr)};
    traccc::detector_conditions_description::data host_det_cond_data{
        vecmem::get_data(host_det_cond)};
    traccc::detector_design_description::buffer device_det_descr{
        [&]() {
            std::vector<unsigned int> sizes(host_det_descr.size());
            for (std::size_t i = 0; i < host_det_descr.size(); ++i) {
                auto this_design = host_det_descr.at(i);

                sizes[i] = std::max(static_cast<unsigned int>(
                                        ((this_design.bin_edges_x()).size())),
                                    static_cast<unsigned int>(
                                        ((this_design.bin_edges_y()).size())));
            }
            return sizes;
        }(),
        device_mr, &host_mr, vecmem::data::buffer_type::resizable};
    copy.setup(device_det_descr)->wait();
    copy(vecmem::get_data(host_det_descr), device_det_descr)->wait();

    traccc::detector_conditions_description::buffer device_det_cond{
        static_cast<traccc::detector_conditions_description::buffer::size_type>(
            host_det_cond.size()),
        device_mr};
    copy.setup(device_det_cond)->ignore();
    copy(vecmem::get_data(host_det_cond), device_det_cond,
         vecmem::copy::type::host_to_device)
        ->ignore();

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector host_det;
    traccc::io::read_detector(host_det, host_mr, detector_opts.detector_file,
                              detector_opts.material_file,
                              detector_opts.grid_file);

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(host_det, device_mr, copy);

    // Output stats
    uint64_t n_cells = 0;
    uint64_t n_measurements = 0;
    uint64_t n_measurements_alpaka = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_spacepoints_alpaka = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_alpaka = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_found_tracks_alpaka = 0;
    uint64_t n_fitted_tracks = 0;
    uint64_t n_fitted_tracks_alpaka = 0;

    // Type definitions
    using host_spacepoint_formation_algorithm =
        traccc::host::silicon_pixel_spacepoint_formation_algorithm;
    using device_spacepoint_formation_algorithm =
        traccc::alpaka::silicon_pixel_spacepoint_formation_algorithm;

    using host_finding_algorithm =
        traccc::host::combinatorial_kalman_filter_algorithm;
    using device_finding_algorithm =
        traccc::alpaka::combinatorial_kalman_filter_algorithm;

    using host_fitting_algorithm = traccc::host::kalman_fitting_algorithm;
    using device_fitting_algorithm = traccc::alpaka::kalman_fitting_algorithm;

    // Algorithm configuration(s).
    const traccc::seedfinder_config seedfinder_config(seeding_opts);
    const traccc::seedfilter_config seedfilter_config(seeding_opts);
    const traccc::spacepoint_grid_config spacepoint_grid_config(seeding_opts);

    detray::propagation::config propagation_config(propagation_opts);

    traccc::finding_config finding_cfg(finding_opts);
    finding_cfg.propagation = propagation_config;

    traccc::fitting_config fitting_cfg(fitting_opts);
    fitting_cfg.propagation = propagation_config;

    // Constant B field for the track finding and fitting
    const traccc::vector3 field_vec(seeding_opts);
    const auto field = traccc::details::make_magnetic_field(bfield_opts);

    traccc::host::clusterization_algorithm ca(
        host_mr, logger().clone("HostClusteringAlg"));
    host_spacepoint_formation_algorithm sf(
        host_mr, logger().clone("HostSpFormationAlg"));
    traccc::host::seeding_algorithm sa(
        seedfinder_config, spacepoint_grid_config, seedfilter_config, host_mr,
        logger().clone("HostSeedingAlg"));
    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(
        track_params_estimation_config, host_mr,
        logger().clone("HostTrackParEstAlg"));
    host_finding_algorithm finding_alg(finding_cfg, host_mr,
                                       logger().clone("HostFindingAlg"));
    host_fitting_algorithm fitting_alg(fitting_cfg, host_mr, host_copy,
                                       logger().clone("HostFittingAlg"));

    traccc::alpaka::clusterization_algorithm ca_alpaka(
        mr, copy, queue, clusterization_opts,
        logger().clone("AlpakaClusteringAlg"));
    traccc::alpaka::measurement_sorting_algorithm ms_alpaka(
        mr, copy, queue, logger().clone("AlpakaMeasSortingAlg"));
    device_spacepoint_formation_algorithm sf_alpaka(
        mr, copy, queue, logger().clone("AlpakaSpFormationAlg"));
    traccc::alpaka::triplet_seeding_algorithm sa_alpaka(
        seedfinder_config, spacepoint_grid_config, seedfilter_config, mr, copy,
        queue, logger().clone("AlpakaSeedingAlg"));
    traccc::alpaka::seed_parameter_estimation_algorithm tp_alpaka(
        track_params_estimation_config, mr, copy, queue,
        logger().clone("AlpakaTrackParEstAlg"));
    device_finding_algorithm finding_alg_alpaka(
        finding_cfg, mr, copy, queue, logger().clone("AlpakaFindingAlg"));
    device_fitting_algorithm fitting_alg_alpaka(
        fitting_cfg, mr, copy, queue, logger().clone("AlpakaFittingAlg"));

    // performance writer
    traccc::seeding_performance_writer sd_performance_writer(
        traccc::seeding_performance_writer::config{},
        logger().clone("SeedingPerformanceWriter"));

    traccc::performance::timing_info elapsedTimes;

    // Loop over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Instantiate host containers/collections
        traccc::host::clusterization_algorithm::output_type
            measurements_per_event{host_mr};
        host_spacepoint_formation_algorithm::output_type spacepoints_per_event{
            host_mr};
        traccc::host::seeding_algorithm::output_type seeds{host_mr};
        traccc::host::track_params_estimation::output_type params{&host_mr};
        host_finding_algorithm::output_type track_candidates{host_mr};
        host_fitting_algorithm::output_type track_states{host_mr};

        // Instantiate alpaka containers/collections
        traccc::edm::measurement_collection::buffer measurements_alpaka_buffer;
        traccc::edm::spacepoint_collection::buffer spacepoints_alpaka_buffer;
        traccc::edm::seed_collection::buffer seeds_alpaka_buffer;
        traccc::bound_track_parameters_collection_types::buffer
            params_alpaka_buffer(0, *mr.host);
        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_candidates_buffer;
        traccc::edm::track_container<traccc::default_algebra>::buffer
            track_states_buffer;

        {
            traccc::performance::timer wall_t("Wall time", elapsedTimes);

            traccc::edm::silicon_cell_collection::host cells_per_event{host_mr};

            {
                traccc::performance::timer t("File reading  (cpu)",
                                             elapsedTimes);
                // Read the cells from the relevant event file into host memory.
                static constexpr bool DEDUPLICATE = true;
                traccc::io::read_cells(
                    cells_per_event, event, input_opts.directory,
                    logger().clone(), &host_det_cond, input_opts.format,
                    DEDUPLICATE, input_opts.use_acts_geom_source);
            }  // stop measuring file reading timer

            n_cells += cells_per_event.size();

            // Create device copy of input collections
            traccc::edm::silicon_cell_collection::buffer cells_buffer(
                static_cast<unsigned int>(cells_per_event.size()), mr.main);
            copy.setup(cells_buffer)->wait();
            copy(vecmem::get_data(cells_per_event), cells_buffer)->wait();

            // Alpaka
            {
                traccc::performance::timer t("Clusterization (alpaka)",
                                             elapsedTimes);
                // Reconstruct it into spacepoints on the device.
                auto unsorted_measurements =
                    ca_alpaka(cells_buffer, device_det_descr, device_det_cond);
                measurements_alpaka_buffer = ms_alpaka(unsorted_measurements);
                queue.synchronize();
            }  // stop measuring clusterization alpaka timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Clusterization  (cpu)",
                                             elapsedTimes);
                measurements_per_event =
                    ca(vecmem::get_data(cells_per_event), host_det_descr_data,
                       host_det_cond_data);
            }  // stop measuring clusterization cpu timer

            // Perform seeding, track finding and fitting only when using a
            // Detray geometry.

            // Alpaka
            {
                traccc::performance::timer t("Spacepoint formation (alpaka)",
                                             elapsedTimes);
                spacepoints_alpaka_buffer =
                    sf_alpaka(detector_buffer, measurements_alpaka_buffer);
                queue.synchronize();
            }  // stop measuring spacepoint formation alpaka timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Spacepoint formation  (cpu)",
                                             elapsedTimes);
                spacepoints_per_event =
                    sf(host_det, vecmem::get_data(measurements_per_event));
            }  // stop measuring spacepoint formation cpu timer

            // Alpaka
            {
                traccc::performance::timer t("Seeding (alpaka)", elapsedTimes);
                seeds_alpaka_buffer = sa_alpaka(spacepoints_alpaka_buffer);
                queue.synchronize();
            }  // stop measuring seeding alpaka timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Seeding  (cpu)", elapsedTimes);
                seeds = sa(vecmem::get_data(spacepoints_per_event));
            }  // stop measuring seeding cpu timer

            // Alpaka
            {
                traccc::performance::timer t("Track params (alpaka)",
                                             elapsedTimes);
                params_alpaka_buffer =
                    tp_alpaka(field, measurements_alpaka_buffer,
                              spacepoints_alpaka_buffer, seeds_alpaka_buffer);
                queue.synchronize();
            }  // stop measuring track params timer

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer t("Track params  (cpu)",
                                             elapsedTimes);
                params = tp(vecmem::get_data(measurements_per_event),
                            vecmem::get_data(spacepoints_per_event),
                            vecmem::get_data(seeds), field_vec);
            }  // stop measuring track params cpu timer

            // Alpaka
            {
                traccc::performance::timer timer{"Track finding (alpaka)",
                                                 elapsedTimes};
                track_candidates_buffer = finding_alg_alpaka(
                    detector_buffer, field, measurements_alpaka_buffer,
                    params_alpaka_buffer);
            }

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer timer{"Track finding (cpu)",
                                                 elapsedTimes};
                track_candidates = finding_alg(
                    host_det, field, vecmem::get_data(measurements_per_event),
                    vecmem::get_data(params));
            }

            // Alpaka
            {
                traccc::performance::timer timer{"Track fitting (alpaka)",
                                                 elapsedTimes};
                track_states_buffer = fitting_alg_alpaka(
                    detector_buffer, field, track_candidates_buffer);
            }

            // CPU
            if (accelerator_opts.compare_with_cpu) {
                traccc::performance::timer timer{"Track fitting (cpu)",
                                                 elapsedTimes};
                track_states = fitting_alg(
                    host_det, field,
                    traccc::edm::track_container<
                        traccc::default_algebra>::const_data(track_candidates));
            }
        }  // Stop measuring wall time

        /*----------------------------------
          compare cpu and alpaka result
          ----------------------------------*/

        traccc::edm::measurement_collection::host measurements_per_event_alpaka{
            host_mr};
        traccc::edm::spacepoint_collection::host spacepoints_per_event_alpaka{
            host_mr};
        traccc::edm::seed_collection::host seeds_alpaka{host_mr};
        traccc::bound_track_parameters_collection_types::host params_alpaka{
            &host_mr};
        traccc::edm::track_collection<traccc::default_algebra>::host
            track_candidates_alpaka{host_mr};
        traccc::edm::track_container<traccc::default_algebra>::host
            track_states_alpaka{host_mr};

        copy(measurements_alpaka_buffer, measurements_per_event_alpaka)->wait();
        copy(spacepoints_alpaka_buffer, spacepoints_per_event_alpaka)->wait();
        copy(seeds_alpaka_buffer, seeds_alpaka)->wait();
        copy(params_alpaka_buffer, params_alpaka)->wait();
        copy(track_candidates_buffer.tracks, track_candidates_alpaka)->wait();
        copy(track_states_buffer.tracks, track_states_alpaka.tracks)->wait();
        copy(track_states_buffer.states, track_states_alpaka.states)->wait();
        queue.synchronize();

        if (accelerator_opts.compare_with_cpu) {

            // Show which event we are currently presenting the results for.
            TRACCC_INFO("===>>> Event " << event << " <<<===");

            // Compare the measurements made on the host and on the device.
            traccc::soa_comparator<traccc::edm::measurement_collection>
                compare_measurements{"measurements"};
            compare_measurements(
                vecmem::get_data(measurements_per_event),
                vecmem::get_data(measurements_per_event_alpaka));

            // Compare the spacepoints made on the host and on the device.
            traccc::soa_comparator<traccc::edm::spacepoint_collection>
                compare_spacepoints{"spacepoints"};
            compare_spacepoints(vecmem::get_data(spacepoints_per_event),
                                vecmem::get_data(spacepoints_per_event_alpaka));

            // Compare the seeds made on the host and on the device
            traccc::soa_comparator<traccc::edm::seed_collection> compare_seeds{
                "seeds", traccc::details::comparator_factory<
                             traccc::edm::seed_collection::const_device::
                                 const_proxy_type>{
                             vecmem::get_data(spacepoints_per_event),
                             vecmem::get_data(spacepoints_per_event_alpaka)}};
            compare_seeds(vecmem::get_data(seeds),
                          vecmem::get_data(seeds_alpaka));

            // Compare the track parameters made on the host and on the device.
            traccc::collection_comparator<traccc::bound_track_parameters<>>
                compare_track_parameters{"track parameters"};
            compare_track_parameters(vecmem::get_data(params),
                                     vecmem::get_data(params_alpaka));

            // Compare tracks found on the host and on the device.
            traccc::soa_comparator<
                traccc::edm::track_collection<traccc::default_algebra>>
                compare_track_candidates{
                    "track candidates",
                    traccc::details::comparator_factory<
                        traccc::edm::track_collection<traccc::default_algebra>::
                            const_device::const_proxy_type>{
                        vecmem::get_data(measurements_per_event),
                        vecmem::get_data(measurements_per_event_alpaka),
                        {},
                        {}}};
            compare_track_candidates(vecmem::get_data(track_candidates.tracks),
                                     vecmem::get_data(track_candidates_alpaka));

            // Compare tracks fitted on the host and on the device.
            traccc::soa_comparator<
                traccc::edm::track_collection<traccc::default_algebra>>
                compare_track_fits{
                    "track fits",
                    traccc::details::comparator_factory<
                        traccc::edm::track_collection<traccc::default_algebra>::
                            const_device::const_proxy_type>{
                        vecmem::get_data(measurements_per_event),
                        vecmem::get_data(measurements_per_event_alpaka),
                        vecmem::get_data(track_states.states),
                        vecmem::get_data(track_states_alpaka.states)}};
            compare_track_fits(vecmem::get_data(track_states.tracks),
                               vecmem::get_data(track_states_alpaka.tracks));
        }
        /// Statistics
        n_measurements += measurements_per_event.size();
        n_spacepoints += spacepoints_per_event.size();
        n_seeds += seeds.size();
        n_measurements_alpaka += measurements_per_event_alpaka.size();
        n_spacepoints_alpaka += spacepoints_per_event_alpaka.size();
        n_seeds_alpaka += seeds_alpaka.size();
        n_found_tracks += track_candidates.tracks.size();
        n_found_tracks_alpaka += track_candidates_alpaka.size();
        n_fitted_tracks += track_states.tracks.size();
        n_fitted_tracks_alpaka += track_states_alpaka.tracks.size();

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
    }

    TRACCC_INFO("==> Statistics ... ");
    TRACCC_INFO("- read    " << n_cells << " cells");
    TRACCC_INFO("- created (cpu)  " << n_measurements << " measurements     ");
    TRACCC_INFO("- created (alpaka)  " << n_measurements_alpaka
                                       << " measurements     ");
    TRACCC_INFO("- created (cpu)  " << n_spacepoints << " spacepoints     ");
    TRACCC_INFO("- created (alpaka) " << n_spacepoints_alpaka
                                      << " spacepoints     ");

    TRACCC_INFO("- created  (cpu) " << n_seeds << " seeds");
    TRACCC_INFO("- created (alpaka) " << n_seeds_alpaka << " seeds");
    TRACCC_INFO("- found (cpu)    " << n_found_tracks << " tracks");
    TRACCC_INFO("- found (alpaka)   " << n_found_tracks_alpaka << " tracks");
    TRACCC_INFO("- fitted (cpu)   " << n_fitted_tracks << " tracks");
    TRACCC_INFO("- fitted (alpaka)  " << n_fitted_tracks_alpaka << " tracks");
    TRACCC_INFO("==>Elapsed times... " << elapsedTimes);

    return 0;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSeqAlpaka", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::clusterization clusterization_opts;
    traccc::opts::track_seeding seeding_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::accelerator accelerator_opts;
    traccc::opts::program_options program_opts{
        "Full Tracking Chain Using Alpaka",
        {detector_opts, bfield_opts, input_opts, clusterization_opts,
         seeding_opts, finding_opts, propagation_opts, performance_opts,
         fitting_opts, accelerator_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(detector_opts, bfield_opts, input_opts, clusterization_opts,
                   seeding_opts, finding_opts, propagation_opts, fitting_opts,
                   performance_opts, accelerator_opts, logger->clone());
}
