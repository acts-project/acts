/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/fitting/kalman_filter/kalman_actor.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/utils/event_data.hpp"
#include "traccc/utils/fill_track_container.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/propagation.hpp"

// Options
#include "traccc/options/detector.hpp"
#include "traccc/options/generation.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_propagation.hpp"

// Performance include(s).
#include "traccc/performance/kalman_filter_comparison.hpp"

// Detray include(s)
#include <detray/io/frontend/detector_reader.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <filesystem>

int main(int argc, char* argv[]) {

    using detector_t = traccc::default_detector::host;
    using algebra_t = typename detector_t::algebra_type;
    using vector3_t = typename detector_t::vector3_type;

    using b_field_t =
        covfie::field<traccc::const_bfield_backend_t<traccc::scalar>>;
    using sf_candidate_t =
        traccc::propagation_validator::candidate_type<detector_t>;

    std::unique_ptr<const traccc::Logger> ilogger = traccc::getDefaultLogger(
        "KalmanFilterValidationCPU", traccc::Logging::Level::INFO);
    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::magnetic_field bfield_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::generation generation_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::program_options program_opts{
        "Kalman Filter validation on the Host",
        {detector_opts, input_opts, generation_opts, propagation_opts},
        argc,
        argv,
        logger().cloneWithSuffix("Options")};

    // Memory resource used by the application.
    vecmem::host_memory_resource host_mr;

    // Truth particle cuts
    constexpr traccc::scalar min_p{50.f * traccc::unit<traccc::scalar>::MeV};
    constexpr traccc::scalar max_r{75.f * traccc::unit<traccc::scalar>::mm};

    TRACCC_INFO("Reading detector from file");

    // Set up the detector reader configuration.
    detray::io::detector_reader_config reader_cfg;
    reader_cfg.add_file(
        traccc::io::get_absolute_path(detector_opts.detector_file));
    if (!detector_opts.material_file.empty()) {
        reader_cfg.add_file(
            traccc::io::get_absolute_path(detector_opts.material_file));
    }
    if (!detector_opts.grid_file.empty()) {
        reader_cfg.add_file(
            traccc::io::get_absolute_path(detector_opts.grid_file));
    }

    // Read the detector.
    auto [io_det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    traccc::host_detector host_det{};
    host_det
        .template set<traccc::detector_traits<typename detector_t::metadata>>(
            std::move(io_det));
    const auto& det = host_det.template as<
        traccc::detector_traits<typename detector_t::metadata>>();

    // Create B-field
    const vector3_t B{0.f, 0.f, 2.f * traccc::unit<traccc::scalar>::T};
    traccc::magnetic_field field = traccc::construct_const_bfield(B);
    // const auto field = traccc::details::make_magnetic_field(bfield_opts);

    TRACCC_INFO("Preparing input data");

    // Check input dir
    std::filesystem::path data_dir{traccc::io::data_directory()};
    std::filesystem::path input_dir = data_dir / input_opts.directory;

    // No existing truth data: Run fast sim
    if (!std::filesystem::exists(input_dir) ||
        std::filesystem::is_empty(input_dir)) {

        TRACCC_INFO("Input data does not exist");
        return EXIT_FAILURE;
    } else {
        TRACCC_INFO("Reading truth data in " << input_opts.directory << "\n");
    }

    detray::propagation_validation_config test_cfg{};
    traccc::seed_generator<detector_t>::config smearing_cfg{};

    // Specific config for the navigation test
    test_cfg.propagation = propagation_opts;

    // Configure the material interaction
    test_cfg.particle = generation_opts.ptc_type;
    test_cfg.do_multiple_scattering = generation_opts.do_multiple_scattering;
    test_cfg.do_energy_loss = generation_opts.do_energy_loss;

    // Prepare the data for the test

    // Initial track parameters from truth particle
    std::vector<traccc::free_track_parameters<algebra_t>> tracks{};
    // The traces of truth hits forward and in reverse order
    std::vector<vecmem::vector<sf_candidate_t>> truth_traces_fw{};
    // Measurements
    typename traccc::edm::measurement_collection::host measurements{host_mr};
    // Collection for bound track parameters and track state/measurement links
    // Input/results of the KF
    traccc::edm::track_container<algebra_t>::host track_container{host_mr};

    traccc::fill_track_containers(ilogger->clone(), &host_det, input_dir,
                                  static_cast<unsigned int>(input_opts.events),
                                  input_opts.use_acts_geom_source, min_p, max_r,
                                  tracks, truth_traces_fw, measurements,
                                  track_container);

    // Save the original truth traces before dummy records are inserted for
    // missing sensitive surfaces
    auto truth_traces_fw_KF = truth_traces_fw;

    std::optional<b_field_t::view_t> field_view{field.template as_field<
        traccc::const_bfield_backend_t<traccc::scalar>>()};

    // Run the test without KF
    detray::propagation_validation(det, names, field_view, test_cfg, tracks,
                                   truth_traces_fw);

    // Set up the containers for the Kalman actor
    traccc::edm::measurement_collection::const_device device_measurements{
        vecmem::get_data(measurements)};

    traccc::kalman_filter_comparison(
        det, names, field, test_cfg, smearing_cfg, ilogger->clone(), tracks,
        truth_traces_fw_KF, device_measurements, track_container);
}
