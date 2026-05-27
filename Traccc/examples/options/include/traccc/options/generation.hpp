/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/options/details/value_array.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/ranges.hpp"

// System include(s).
#include <cstddef>

namespace traccc::opts {

/// Configuration for particle / event generation
class generation : public interface {

    public:
    /// @name Configurable options
    /// @{

    /// The number of events to generate
    std::size_t events = 1;

    /// The number of particles to generate per event
    unsigned int gen_nparticles{1u};
    /// Vertex position [mm]
    opts::value_array<float, 3> vertex{0.f, 0.f, 0.f};
    /// Standard deviation of the vertex position [mm]
    opts::value_array<float, 3> vertex_stddev{0.f, 0.f, 0.f};
    /// Range of momentum [GeV]
    opts::value_array<float, 2> mom_range{1.f, 1.f};
    /// Range of phi [rad]
    opts::value_array<float, 2> phi_range{-180.f, 180.f};
    /// Range of eta
    opts::value_array<float, 2> eta_range{-2.f, 2.f};
    /// Range of theta [rad] (corresponding to eta range of [-2,2])
    opts::value_array<float, 2> theta_range{eta_to_theta_range(eta_range)};
    /// Whether to include multiple scattering
    bool do_multiple_scattering{true};
    /// Whether to include energy loss
    bool do_energy_loss{true};

    /// PDG number for particle type (Default: muon)
    int pdg_number = 13;

    /// @}

    /// @name Derived options
    /// @{

    /// Particle type
    traccc::pdg_particle<traccc::scalar> ptc_type =
        traccc::muon<traccc::scalar>();

    /// @}

    /// Constructor
    generation();

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct generation

}  // namespace traccc::opts
