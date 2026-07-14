/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/generation.hpp"
#include "traccc/options/output_data.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/measurement_smearer.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/simulation/smearing_writer.hpp"
#include "traccc/utils/propagation.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Boost include(s).
#include <boost/filesystem.hpp>

using namespace traccc;

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSimulate", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::detector det_opts;
    traccc::opts::generation generation_opts;
    traccc::opts::output_data output_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::program_options program_opts{
        "Detector Simulation",
        {det_opts, generation_opts, output_opts, propagation_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    /// Type declarations
    using host_detector_type = traccc::default_detector::host;
    using uniform_gen_t =
        detray::detail::random_numbers<scalar,
                                       std::uniform_real_distribution<scalar>>;
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters<>,
                                       uniform_gen_t>;

    // B field value and its type
    // @TODO: Set B field as argument
    using b_field_t =
        covfie::field<traccc::const_bfield_backend_t<traccc::scalar>>;
    const traccc::vector3 B{0, 0, 2 * traccc::unit<traccc::scalar>::T};
    const b_field_t field =
        traccc::construct_const_bfield(B)
            .as_field<traccc::const_bfield_backend_t<traccc::scalar>>();

    // Read the detector
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(traccc::io::get_absolute_path(det_opts.detector_file));
    if (!det_opts.material_file.empty()) {
        reader_cfg.add_file(
            traccc::io::get_absolute_path(det_opts.material_file));
    }
    if (!det_opts.grid_file.empty()) {
        reader_cfg.add_file(traccc::io::get_absolute_path(det_opts.grid_file));
    }

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    /***************************
     * Generate simulation data
     ***************************/

    // Origin of particles
    generator_type::configuration gen_cfg{};
    gen_cfg.n_tracks(generation_opts.gen_nparticles);
    gen_cfg.origin(traccc::point3{generation_opts.vertex[0],
                                  generation_opts.vertex[1],
                                  generation_opts.vertex[2]});
    gen_cfg.origin_stddev(traccc::point3{generation_opts.vertex_stddev[0],
                                         generation_opts.vertex_stddev[1],
                                         generation_opts.vertex_stddev[2]});
    gen_cfg.phi_range(generation_opts.phi_range);
    gen_cfg.theta_range(generation_opts.theta_range);
    gen_cfg.mom_range(generation_opts.mom_range);
    generator_type generator(gen_cfg);

    // Smearing value for measurements
    traccc::measurement_smearer<traccc::default_algebra> meas_smearer(
        50 * traccc::unit<scalar>::um, 50 * traccc::unit<scalar>::um);

    using writer_type = traccc::smearing_writer<
        traccc::measurement_smearer<traccc::default_algebra>>;

    // Writer config
    typename writer_type::config smearer_writer_cfg{meas_smearer};

    // Run simulator
    const std::string full_path = io::data_directory() + output_opts.directory;

    boost::filesystem::create_directories(full_path);

    auto sim = traccc::simulator<host_detector_type, b_field_t, generator_type,
                                 writer_type>(
        generation_opts.ptc_type, generation_opts.events, host_det, field,
        std::move(generator), std::move(smearer_writer_cfg), full_path);

    sim.get_config().propagation = propagation_opts;

    sim.run();

    return 1;
}
