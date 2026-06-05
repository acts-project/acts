/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/generation.hpp"
#include "traccc/options/output_data.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/telescope_detector.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/measurement_smearer.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/simulation/smearing_writer.hpp"

// detray include(s).
#include <detray/geometry/mask.hpp>
#include <detray/geometry/shapes/rectangle2D.hpp>
#include <detray/material/material.hpp>
#include <detray/test/common/build_telescope_detector.hpp>
#include <detray/tracks/ray.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Boost include(s).
#include <boost/filesystem.hpp>

using namespace traccc;

int simulate(const traccc::opts::generation& generation_opts,
             const traccc::opts::output_data& output_opts,
             const traccc::opts::track_propagation& propagation_opts,
             const traccc::opts::telescope_detector& telescope_opts) {

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        detray::detail::random_numbers<scalar,
                                       std::uniform_real_distribution<scalar>>;

    // Memory resource
    vecmem::host_memory_resource host_mr;

    /*****************************
     * Build a telescope geometry
     *****************************/

    const auto align_vec = telescope_opts.align_vector;
    vector3 align_axis{align_vec[0u], align_vec[1u], align_vec[2u]};
    align_axis = vector::normalize(align_axis);

    // Plane alignment direction (aligned to z-axis)
    detray::detail::ray<traccc::default_algebra> traj{
        {0, 0, 0}, 0, align_axis, -1};
    // Position of planes (in mm unit)
    std::vector<scalar> plane_positions;

    for (unsigned int i = 0; i < telescope_opts.n_planes; i++) {
        plane_positions.push_back(static_cast<float>(i) *
                                  telescope_opts.spacing);
    }

    // B field value and its type
    using b_field_t = covfie::field<traccc::const_bfield_backend_t<scalar>>;
    const vector3 B{0, 0, 2 * traccc::unit<scalar>::T};
    b_field_t field = traccc::construct_const_bfield(B)
                          .as_field<traccc::const_bfield_backend_t<scalar>>();

    // Set material and thickness
    detray::material<scalar> mat;
    if (!telescope_opts.empty_material) {
        mat = detray::silicon<scalar>();
    } else {
        mat = detray::vacuum<scalar>();
    }

    const scalar thickness = telescope_opts.thickness;

    // Use rectangle surfaces
    detray::mask<detray::rectangle2D, traccc::default_algebra> rectangle{
        0u, telescope_opts.half_length, telescope_opts.half_length};

    detray::tel_det_config tel_cfg{rectangle};
    tel_cfg.positions(plane_positions);
    tel_cfg.module_material(mat);
    tel_cfg.mat_thickness(thickness);
    tel_cfg.pilot_track(traj);

    const auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    /***************************
     * Generate simulation data
     ***************************/

    // Origin of particles
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters<>,
                                       uniform_gen_t>;
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
    gen_cfg.charge(generation_opts.ptc_type.charge());
    generator_type generator(gen_cfg);

    // Smearing value for measurements
    traccc::measurement_smearer<traccc::default_algebra> meas_smearer(
        telescope_opts.smearing, telescope_opts.smearing);

    // Type declarations
    using detector_type = decltype(det);
    using writer_type = traccc::smearing_writer<
        traccc::measurement_smearer<traccc::default_algebra>>;

    // Writer config
    typename writer_type::config smearer_writer_cfg{meas_smearer};

    // Run simulator
    const std::string full_path = io::data_directory() + output_opts.directory;

    boost::filesystem::create_directories(full_path);

    auto sim = traccc::simulator<detector_type, b_field_t, generator_type,
                                 writer_type>(
        generation_opts.ptc_type, generation_opts.events, det, field,
        std::move(generator), std::move(smearer_writer_cfg), full_path);
    sim.get_config().propagation = propagation_opts;

    sim.run();

    // Create detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, name_map, writer_cfg);

    return 1;
}

// The main routine
//
int main(int argc, char* argv[]) {
    std::unique_ptr<const traccc::Logger> logger = traccc::getDefaultLogger(
        "TracccExampleSimulateTelescope", traccc::Logging::Level::INFO);

    // Program options.
    traccc::opts::generation generation_opts;
    traccc::opts::output_data output_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::telescope_detector telescope_opts;
    traccc::opts::program_options program_opts{
        "Telescope-Detector Simulation",
        {generation_opts, output_opts, propagation_opts, telescope_opts},
        argc,
        argv,
        logger->cloneWithSuffix("Options")};

    // Run the application.
    return simulate(generation_opts, output_opts, propagation_opts,
                    telescope_opts);
}
