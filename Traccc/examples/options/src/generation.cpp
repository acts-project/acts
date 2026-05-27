/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/generation.hpp"

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/ranges.hpp"

// System include(s).
#include <sstream>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

generation::generation() : interface("Particle Generation Options") {

    m_desc.add_options()("gen-events",
                         po::value(&events)->default_value(events),
                         "The number of events to generate");
    m_desc.add_options()(
        "gen-nparticles",
        po::value(&gen_nparticles)->default_value(gen_nparticles),
        "The number of particles to generate per event");
    m_desc.add_options()(
        "gen-vertex-xyz-mm",
        po::value(&vertex)->value_name("X:Y:Z")->default_value(vertex),
        "Vertex [mm]");
    m_desc.add_options()("gen-vertex-xyz-std-mm",
                         po::value(&vertex_stddev)
                             ->value_name("X:Y:Z")
                             ->default_value(vertex_stddev),
                         "Standard deviation of the vertex [mm]");
    m_desc.add_options()(
        "gen-mom-gev",
        po::value(&mom_range)->value_name("MIN:MAX")->default_value(mom_range),
        "Range of momentum [GeV]");
    m_desc.add_options()(
        "gen-phi-degree",
        po::value(&phi_range)->value_name("MIN:MAX")->default_value(phi_range),
        "Range of phi [Degree]");
    m_desc.add_options()(
        "gen-eta",
        po::value(&eta_range)->value_name("MIN:MAX")->default_value(eta_range),
        "Range of eta");
    m_desc.add_options()("gen-theta",
                         po::value(&theta_range)
                             ->value_name("MIN:MAX")
                             ->default_value(theta_range),
                         "Range of theta in degree");
    m_desc.add_options()("particle-type",
                         po::value(&pdg_number)->default_value(pdg_number),
                         "PDG number for the particle type");
}

void generation::read(const po::variables_map &vm) {

    vertex *= traccc::unit<float>::mm;
    vertex_stddev *= traccc::unit<float>::mm;
    mom_range *= traccc::unit<float>::GeV;
    phi_range *= traccc::unit<float>::degree;

    // The eta and theta range can not be specified at the same time
    if (vm.count("gen-eta") && !vm["gen-eta"].defaulted() &&
        vm.count("gen-theta") && !vm["gen-theta"].defaulted()) {
        throw std::logic_error(
            std::string("Conflicting options 'gen-eta' and 'gen-theta'"));
    } else if (vm.count("gen-eta") && !vm["gen-eta"].defaulted()) {
        theta_range = eta_to_theta_range(eta_range);
    } else if (vm.count("gen-theta") && !vm["gen-theta"].defaulted()) {
        theta_range *= traccc::unit<float>::degree;
        eta_range = theta_to_eta_range(theta_range);
    }

    ptc_type = detail::particle_from_pdg_number<traccc::scalar>(pdg_number);
}

std::unique_ptr<configuration_printable> generation::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of events", std::to_string(events)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of particles", std::to_string(gen_nparticles)));
    std::ostringstream vertex_ss;
    vertex_ss << vertex / traccc::unit<float>::mm << " mm";
    cat->add_child(
        std::make_unique<configuration_kv_pair>("Vertex", vertex_ss.str()));
    std::ostringstream vertex_dev_ss;
    vertex_dev_ss << vertex_stddev / traccc::unit<float>::mm << " mm";
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Vertex standard deviation", vertex_dev_ss.str()));
    std::ostringstream mom_range_ss;
    mom_range_ss << mom_range / traccc::unit<float>::GeV << " GeV";
    cat->add_child(std::make_unique<configuration_kv_pair>("Momentum range",
                                                           mom_range_ss.str()));
    std::ostringstream phi_range_ss;
    phi_range_ss << phi_range / traccc::unit<float>::degree << " deg";
    cat->add_child(std::make_unique<configuration_kv_pair>("Phi range",
                                                           phi_range_ss.str()));
    std::ostringstream eta_range_ss;
    eta_range_ss << eta_range;
    cat->add_child(std::make_unique<configuration_kv_pair>("Eta range",
                                                           eta_range_ss.str()));
    std::ostringstream theta_range_ss;
    theta_range_ss << theta_range / traccc::unit<float>::degree << " deg";
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Theta range", theta_range_ss.str()));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "PGD number", std::to_string(pdg_number)));

    return cat;
}

}  // namespace traccc::opts
