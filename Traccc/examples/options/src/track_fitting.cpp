/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/track_fitting.hpp"

#include "traccc/examples/utils/printable.hpp"
#include "traccc/utils/particle.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

track_fitting::track_fitting() : interface("Track Fitting Options") {

    m_desc.add_options()(
        "fit-num-iterations",
        po::value(&m_config.n_iterations)->default_value(m_config.n_iterations),
        "Number of iterations for the track fit");
    m_desc.add_options()(
        "fit-particle-hypothesis",
        po::value(&m_pdg)->value_name("PDG")->default_value(m_pdg),
        "Particle hypothesis for the track fit");
    m_desc.add_options()(
        "fit-covariance-inflation-factor",
        po::value(&m_config.covariance_inflation_factor)
            ->default_value(m_config.covariance_inflation_factor),
        "Covariance inflation factor for the track fit");
    m_desc.add_options()(
        "geo_id-sequence-size-factor",
        po::value(&m_config.surface_sequence_size_factor)
            ->default_value(m_config.surface_sequence_size_factor),
        "Size factor for the geometry identifier sequence used in the backward "
        "filter");
    m_desc.add_options()(
        "min-geo_id-sequence-capacity",
        po::value(&m_config.min_surface_sequence_capacity)
            ->default_value(m_config.min_surface_sequence_capacity),
        "Minimum capacity of geometry identifier sequence");
    m_desc.add_options()(
        "min-total-momentum",
        po::value(&m_config.min_p)->default_value(m_config.min_p),
        "Minimum total track momentum [GeV]");
    m_desc.add_options()(
        "min-transverse-momentum",
        po::value(&m_config.min_pT)->default_value(m_config.min_pT),
        "Minimum transverse track momentum [GeV]");
}

void track_fitting::read(const po::variables_map &) {
    m_config.min_p *= traccc::unit<float>::GeV;
    m_config.min_pT *= traccc::unit<float>::GeV;
}

track_fitting::operator fitting_config() const {

    fitting_config out = m_config;

    out.ptc_hypothesis =
        detail::particle_from_pdg_number<traccc::scalar>(m_pdg);

    return out;
}

std::unique_ptr<configuration_printable> track_fitting::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Number of iterations", std::to_string(m_config.n_iterations)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Particle hypothesis PDG", std::to_string(m_pdg)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Covariance inflation factor",
        std::to_string(m_config.covariance_inflation_factor)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Geometry identifier sequence size factor",
        std::to_string(m_config.surface_sequence_size_factor)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum capacity of geometry identifier sequence",
        std::to_string(m_config.min_surface_sequence_capacity)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum pT",
        std::format("{} GeV", m_config.min_pT / traccc::unit<float>::GeV)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum p",
        std::format("{} GeV", m_config.min_p / traccc::unit<float>::GeV)));

    return cat;
}
}  // namespace traccc::opts
