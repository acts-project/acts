/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/efficiency/finding_performance_writer.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/examples/print_fitted_tracks_statistics.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
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
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/utils/propagation.hpp"
#include "traccc/utils/seed_generator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <cassert>
#include <cstdlib>
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
            const traccc::opts::truth_finding& truth_finding_opts,
            const traccc::opts::track_matching& track_matching_opts,
            std::unique_ptr<const traccc::Logger> ilogger) {
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    // Copy obejct
    vecmem::copy copy;

    // Performance writer
    traccc::finding_performance_writer find_performance_writer(
        traccc::finding_performance_writer::config{
            .truth_config = truth_finding_opts,
            .track_truth_config = track_matching_opts},
        logger().clone("FindingPerformanceWriter"));
    traccc::fitting_performance_writer fit_performance_writer(
        traccc::fitting_performance_writer::config{},
        logger().clone("FittingPerformanceWriter"));

    /*****************************
     * Build a geometry
     *****************************/

    // B field value
    const auto field = traccc::details::make_magnetic_field(bfield_opts);

    // Construct a Detray detector object, if supported by the configuration.
    traccc::host_detector polymorphic_detector;
    traccc::io::read_detector(
        polymorphic_detector, host_mr, detector_opts.detector_file,
        detector_opts.material_file, detector_opts.grid_file);

    /*****************************
     * Do the reconstruction
     *****************************/

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
    typename traccc::finding_config cfg(finding_opts);
    cfg.propagation = propagation_config;

    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(
        cfg, host_mr, logger().clone("FindingAlg"));

    // Fitting algorithm object
    traccc::fitting_config fit_cfg(fitting_opts);
    fit_cfg.propagation = propagation_config;

    traccc::host::kalman_fitting_algorithm host_fitting(
        fit_cfg, host_mr, copy, logger().clone("FittingAlg"));

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

        // Seed generator
        host_detector_visitor<detector_type_list>(
            polymorphic_detector,
            [&]<typename detector_traits_t>(
                const typename detector_traits_t::host& det) {
                traccc::seed_generator<typename detector_traits_t::host> sg(
                    det, stddevs);
                evt_data.generate_truth_candidates(
                    truth_track_candidates, truth_measurements, sg, host_mr,
                    truth_finding_opts.m_pT_min);
            });
        truth_track_candidates.measurements =
            vecmem::get_data(truth_measurements);

        // Prepare truth seeds
        traccc::bound_track_parameters_collection_types::host seeds(&host_mr);
        const std::size_t n_tracks = truth_track_candidates.tracks.size();
        for (std::size_t i_trk = 0; i_trk < n_tracks; i_trk++) {
            seeds.push_back(truth_track_candidates.tracks.at(i_trk).params());
        }

        std::cout << "Number of seeds: " << seeds.size() << std::endl;

        // Read measurements
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::io::read_measurements(
            measurements_per_event, event, input_opts.directory,
            (input_opts.use_acts_geom_source ? &polymorphic_detector : nullptr),
            nullptr, nullptr, input_opts.format);

        // Run finding
        auto track_candidates = host_finding(
            polymorphic_detector, field,
            vecmem::get_data(measurements_per_event), vecmem::get_data(seeds));

        std::cout << "Number of found tracks: "
                  << track_candidates.tracks.size() << std::endl;

        // Run fitting
        auto track_states = host_fitting(
            polymorphic_detector, field,
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                track_candidates));

        details::print_fitted_tracks_statistics(track_states, logger());

        const std::size_t n_fitted_tracks = track_states.tracks.size();

        if (performance_opts.run) {
            find_performance_writer.write(
                traccc::edm::track_container<
                    traccc::default_algebra>::const_data(track_candidates),
                evt_data);

            for (std::size_t i = 0; i < n_fitted_tracks; i++) {
                host_detector_visitor<detector_type_list>(
                    polymorphic_detector,
                    [&]<typename detector_traits_t>(
                        const typename detector_traits_t::host& det) {
                        fit_performance_writer.write(
                            track_states.tracks.at(i), track_states.states,
                            measurements_per_event, det, evt_data);
                    });
            }
        }
    }

    if (performance_opts.run) {
        find_performance_writer.finalize();
        fit_performance_writer.finalize();
    }

    return EXIT_SUCCESS;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleTruthFinding", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_finding finding_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::track_fitting fitting_opts;
    traccc::opts::performance performance_opts;
    traccc::opts::truth_finding truth_finding_config;
    traccc::opts::track_matching track_matching_opts;
    traccc::opts::program_options program_opts{
        "Truth Track Finding on the Host",
        {detector_opts, bfield_opts, input_opts, finding_opts, propagation_opts,
         fitting_opts, performance_opts, truth_finding_config,
         track_matching_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return seq_run(finding_opts, propagation_opts, fitting_opts, input_opts,
                   detector_opts, bfield_opts, performance_opts,
                   truth_finding_config, track_matching_opts, logger->clone());
}
