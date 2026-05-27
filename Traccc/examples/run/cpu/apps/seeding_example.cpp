/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/examples/print_fitted_tracks_statistics.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/propagation.hpp"

// io
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/utils.hpp"

// algorithms
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// performance
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/efficiency/seeding_performance_writer.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"

// options
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
#include "traccc/options/track_resolution.hpp"
#include "traccc/options/track_seeding.hpp"
#include "traccc/options/truth_finding.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace traccc;

int seq_run(const traccc::opts::track_seeding& seeding_opts,
            const traccc::opts::track_finding& finding_opts,
            const traccc::opts::track_propagation& propagation_opts,
            const traccc::opts::track_resolution& resolution_opts,
            const traccc::opts::track_fitting& fitting_opts,
            const traccc::opts::input_data& input_opts,
            const traccc::opts::detector& detector_opts,
            const traccc::opts::magnetic_field& bfield_opts,
            const traccc::opts::performance& performance_opts,
            const traccc::opts::truth_finding& truth_finding_opts,
            const traccc::opts::seed_matching& seed_matching_opts,
            const traccc::opts::track_matching& track_matching_opts,
            std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Copy obejct
    vecmem::copy copy;

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
    traccc::finding_performance_writer::config ar_writer_cfg;
    ar_writer_cfg.file_path = "performance_track_ambiguity_resolution.root";
    ar_writer_cfg.algorithm_name = "ambiguity_resolution";
    traccc::finding_performance_writer ar_performance_writer(
        ar_writer_cfg, logger().clone("AmbiResFindingPerformanceWriter"));
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{},
        logger().clone("FittingPerformanceWriter"));

    // Output stats
    uint64_t n_spacepoints = 0;
    uint64_t n_measurements = 0;
    uint64_t n_seeds = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_ambiguity_free_tracks = 0;
    uint64_t n_fitted_tracks = 0;

    /*****************************
     * Build a geometry
     *****************************/

    // B field value
    const auto field = traccc::details::make_magnetic_field(bfield_opts);
    const traccc::vector3 field_vec(seeding_opts);

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector detector;
    traccc::io::read_detector(detector, host_mr, detector_opts.detector_file,
                              detector_opts.material_file,
                              detector_opts.grid_file);

    // Seeding algorithm
    const traccc::seedfinder_config seedfinder_config(seeding_opts);
    const traccc::seedfilter_config seedfilter_config(seeding_opts);
    const traccc::spacepoint_grid_config spacepoint_grid_config(seeding_opts);
    traccc::host::seeding_algorithm sa(
        seedfinder_config, spacepoint_grid_config, seedfilter_config, host_mr,
        logger().clone("SeedingAlg"));
    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(track_params_estimation_config,
                                             host_mr,
                                             logger().clone("TrackParEstAlg"));

    // Propagation configuration
    detray::propagation::config propagation_config(propagation_opts);

    // Finding algorithm configuration
    traccc::finding_config cfg(finding_opts);
    cfg.propagation = propagation_config;

    traccc::host::combinatorial_kalman_filter_algorithm host_finding(
        cfg, host_mr, logger().clone("FindingAlg"));

    traccc::host::greedy_ambiguity_resolution_algorithm::config_type
        host_ambiguity_config(resolution_opts);
    traccc::host::greedy_ambiguity_resolution_algorithm
        host_ambiguity_resolution(host_ambiguity_config, host_mr,
                                  logger().clone("AmbiguityResolution"));

    // Fitting algorithm object
    traccc::fitting_config fit_cfg(fitting_opts);
    fit_cfg.propagation = propagation_config;

    traccc::host::kalman_fitting_algorithm host_fitting(
        fit_cfg, host_mr, copy, logger().clone("FittingAlg"));

    // Loop over events
    for (std::size_t event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event) {

        // Read the hits from the relevant event file
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::edm::spacepoint_collection::host spacepoints_per_event{host_mr};
        traccc::io::read_spacepoints(
            spacepoints_per_event, measurements_per_event, event,
            input_opts.directory,
            (input_opts.use_acts_geom_source ? &detector : nullptr), nullptr,
            nullptr, input_opts.format);
        n_measurements += measurements_per_event.size();
        n_spacepoints += spacepoints_per_event.size();

        /*----------------
             Seeding
          ---------------*/

        auto seeds = sa(vecmem::get_data(spacepoints_per_event));

        /*----------------------------
           Track Parameter Estimation
          ----------------------------*/

        auto params = tp(vecmem::get_data(measurements_per_event),
                         vecmem::get_data(spacepoints_per_event),
                         vecmem::get_data(seeds), field_vec);

        /*------------------------
           Track Finding with CKF
          ------------------------*/

        auto track_candidates = host_finding(
            detector, field, vecmem::get_data(measurements_per_event),
            vecmem::get_data(params));
        n_found_tracks += track_candidates.tracks.size();

        /*-----------------------------------------
           Ambiguity Resolution with Greedy Solver
          -----------------------------------------*/

        // TODO: Fix me
        /*auto track_candidates_ar = host_ambiguity_resolution(
            traccc::edm::track_container<default_algebra>::const_data(
                track_candidates));
        n_ambiguity_free_tracks += track_candidates_ar.tracks.size();*/

        /*------------------------
           Track Fitting with KF
          ------------------------*/

        auto track_states = host_fitting(
            detector, field,
            traccc::edm::track_container<default_algebra>::const_data(
                track_candidates));
        n_fitted_tracks += track_states.tracks.size();

        /*------------
           Statistics
          ------------*/

        details::print_fitted_tracks_statistics(track_states, logger());
        n_spacepoints += spacepoints_per_event.size();
        n_seeds += seeds.size();

        /*------------
          Writer
          ------------*/

        if (performance_opts.run) {

            traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                        input_opts.use_acts_geom_source,
                                        &detector, input_opts.format, false);

            sd_performance_writer.write(
                vecmem::get_data(seeds),
                vecmem::get_data(spacepoints_per_event),
                vecmem::get_data(measurements_per_event), evt_data);

            find_performance_writer.write(
                traccc::edm::track_container<
                    traccc::default_algebra>::const_data(track_candidates),
                evt_data);

            /*ar_performance_writer.write(
                traccc::edm::track_container<
                    traccc::default_algebra>::const_data(track_candidates_ar),
                evt_data);*/

            for (unsigned int i = 0; i < track_states.tracks.size(); i++) {
                host_detector_visitor<detector_type_list>(
                    detector, [&]<typename detector_traits_t>(
                                  const typename detector_traits_t::host& det) {
                        fit_performance_writer.write(
                            track_states.tracks.at(i), track_states.states,
                            measurements_per_event, det, evt_data);
                    });
            }
        }
    }

    if (performance_opts.run) {
        sd_performance_writer.finalize();
        find_performance_writer.finalize();
        ar_performance_writer.finalize();
        fit_performance_writer.finalize();
    }

    TRACCC_INFO("==> Statistics ... ");
    TRACCC_INFO("- read    " << n_spacepoints << " spacepoints");
    TRACCC_INFO("- read    " << n_measurements << " measurements");
    TRACCC_INFO("- created (cpu)  " << n_seeds << " seeds");
    TRACCC_INFO("- created (cpu)  " << n_found_tracks << " found tracks");
    TRACCC_INFO("- created (cpu)  " << n_ambiguity_free_tracks
                                    << " ambiguity free tracks");
    TRACCC_INFO("- created (cpu)  " << n_fitted_tracks << " fitted tracks");

    return EXIT_SUCCESS;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSeeding", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_seeding seeding_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_resolution resolution_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::truth_finding truth_finding_opts;
    traccc::opts::seed_matching seed_matching_opts;
    traccc::opts::track_matching track_matching_opts;
    traccc::opts::program_options program_opts{
        "Full Tracking Chain on the Host (without clusterization)",
        {detector_opts, bfield_opts, input_opts, seeding_opts, finding_opts,
         propagation_opts, resolution_opts, fitting_opts, performance_opts,
         truth_finding_opts, seed_matching_opts, track_matching_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(seeding_opts, finding_opts, propagation_opts,
                   resolution_opts, fitting_opts, input_opts, detector_opts,
                   bfield_opts, performance_opts, truth_finding_opts,
                   seed_matching_opts, track_matching_opts, logger->clone());
}
