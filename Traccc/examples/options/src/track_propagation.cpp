/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/track_propagation.hpp"

#include "traccc/definitions/common.hpp"
#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>
#include <limits>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

track_propagation::track_propagation()
    : interface("Track Propagation Options") {
    m_search_window[0] = m_config.navigation.search_window[0];
    m_search_window[1] = m_config.navigation.search_window[1];

    m_desc.add_options()("constraint-step-size-mm",
                         po::value(&(m_config.stepping.step_constraint))
                             ->default_value(m_config.stepping.step_constraint /
                                             traccc::unit<float>::mm),
                         "The constrained step size [mm]");
    m_desc.add_options()(
        "overstep-tolerance-um",
        po::value(&(m_config.navigation.intersection.overstep_tolerance))
            ->default_value(
                m_config.navigation.intersection.overstep_tolerance /
                traccc::unit<float>::um),
        "The overstep tolerance [um]");
    m_desc.add_options()(
        "min-mask-tolerance-mm",
        po::value(&(m_config.navigation.intersection.min_mask_tolerance))
            ->default_value(
                m_config.navigation.intersection.min_mask_tolerance /
                traccc::unit<float>::mm),
        "The minimum mask tolerance [mm]");
    m_desc.add_options()(
        "max-mask-tolerance-mm",
        po::value(&(m_config.navigation.intersection.max_mask_tolerance))
            ->default_value(
                m_config.navigation.intersection.max_mask_tolerance /
                traccc::unit<float>::mm),
        "The maximum mask tolerance [mm]");
    m_desc.add_options()(
        "search-window",
        po::value(&m_search_window)->default_value(m_search_window),
        "Size of the grid surface search window");
    m_desc.add_options()(
        "mask-tolerance-scaling",
        po::value(&(m_config.navigation.intersection.mask_tolerance_scalor))
            ->default_value(
                m_config.navigation.intersection.mask_tolerance_scalor),
        "Scale factor between min. and max. mask tolerance with surface "
        "distance");

    m_desc.add_options()(
        "accumulated-noise-factor",
        po::value(&(m_config.navigation.accumulated_error))
            ->default_value(m_config.navigation.accumulated_error),
        "Scale factor on the total track path length to model accumualted "
        "noise [%]");

    m_desc.add_options()(
        "scattering-stddevs",
        po::value(&(m_config.navigation.n_scattering_stddev))
            ->default_value(m_config.navigation.n_scattering_stddev),
        "Number of angle standard deviations to use for the noise modelling");
    m_desc.add_options()("rk-tolerance-mm",
                         po::value(&(m_config.stepping.rk_error_tol))
                             ->default_value(m_config.stepping.rk_error_tol /
                                             traccc::unit<float>::mm),
                         "The Runge-Kutta stepper tolerance [mm]");
    m_desc.add_options()("stepping-min-stepsize",
                         po::value(&(m_config.stepping.min_stepsize))
                             ->default_value(m_config.stepping.min_stepsize /
                                             traccc::unit<float>::mm),
                         "The minimum step size [mm]");
    m_desc.add_options()("stepping-path-limit",
                         po::value(&(m_config.stepping.path_limit))
                             ->default_value(m_config.stepping.path_limit /
                                             traccc::unit<float>::m),
                         "The maximum path length for the stepper [m]");
    m_desc.add_options()("stepping-max-rk-updates",
                         po::value(&(m_config.stepping.max_rk_updates))
                             ->default_value(m_config.stepping.max_rk_updates),
                         "The maximum number of Runge-Kutta updates");
    m_desc.add_options()("stepping-use-mean-loss",
                         po::value(&(m_config.stepping.use_mean_loss))
                             ->default_value(m_config.stepping.use_mean_loss),
                         "Enable the Bethe energy loss model");

    m_desc.add_options()(
        "stepping-do-covariance-transport",
        po::value(&(m_config.stepping.do_covariance_transport))
            ->default_value(m_config.stepping.do_covariance_transport),
        "Enable covariance transport in the stepper");
    m_desc.add_options()(
        "stepping-use-eloss-gradient",
        po::value(&(m_config.stepping.use_eloss_gradient))
            ->default_value(m_config.stepping.use_eloss_gradient),
        "Enable the energy loss gradient in covariance transport");
    m_desc.add_options()(
        "stepping-use-field-gradient",
        po::value(&(m_config.stepping.use_field_gradient))
            ->default_value(m_config.stepping.use_field_gradient),
        "Enable the B-field gradient in covariance transport");
}

void track_propagation::read(const po::variables_map &) {

    m_config.stepping.step_constraint *= traccc::unit<float>::mm;
    m_config.stepping.rk_error_tol *= traccc::unit<float>::mm;
    m_config.navigation.intersection.overstep_tolerance *=
        traccc::unit<float>::um;
    m_config.navigation.intersection.min_mask_tolerance *=
        traccc::unit<float>::mm;
    m_config.navigation.intersection.max_mask_tolerance *=
        traccc::unit<float>::mm;
    m_config.navigation.search_window = m_search_window;
    m_config.navigation.accumulated_error /= 100.f;

    m_config.stepping.min_stepsize *= traccc::unit<float>::mm;
    m_config.stepping.path_limit *= traccc::unit<float>::m;
}

track_propagation::operator detray::propagation::config() const {
    return m_config;
}

std::unique_ptr<configuration_printable> track_propagation::as_printable()
    const {
    auto cat_nav = std::make_unique<configuration_category>("Navigation");

    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Min mask tolerance",
        std::to_string(m_config.navigation.intersection.min_mask_tolerance /
                       traccc::unit<float>::mm) +
            " mm"));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Max mask tolerance",
        std::to_string(m_config.navigation.intersection.max_mask_tolerance /
                       traccc::unit<float>::mm) +
            " mm"));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Mask tolerance scaling",
        std::to_string(
            m_config.navigation.intersection.mask_tolerance_scalor)));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Path tolerance",
        std::to_string(m_config.navigation.intersection.path_tolerance /
                       traccc::unit<float>::um) +
            " um"));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Overstep tolerance",
        std::to_string(m_config.navigation.intersection.overstep_tolerance /
                       traccc::unit<float>::um) +
            " um"));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Search window",
        std::to_string(m_config.navigation.search_window[0]) + " x " +
            std::to_string(m_config.navigation.search_window[1])));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "Scale factor for accumulated noise",
        std::to_string(m_config.navigation.accumulated_error * 100.f) + " %"));
    cat_nav->add_child(std::make_unique<configuration_kv_pair>(
        "# scattering stddevs to assume",
        std::to_string(m_config.navigation.n_scattering_stddev)));

    auto cat_tsp = std::make_unique<configuration_category>("Transport");

    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Min step size", std::to_string(m_config.stepping.min_stepsize /
                                        traccc::unit<float>::mm) +
                             " mm"));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Runge-Kutta tolerance", std::to_string(m_config.stepping.rk_error_tol /
                                                traccc::unit<float>::mm) +
                                     " mm"));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Max step updates", std::to_string(m_config.stepping.max_rk_updates)));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Step size constraint",
        std::to_string(m_config.stepping.step_constraint /
                       traccc::unit<float>::mm) +
            " mm"));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Path limit",
        std::to_string(m_config.stepping.path_limit / traccc::unit<float>::m) +
            " m"));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Enable Bethe energy loss",
        std::format("{}", m_config.stepping.use_mean_loss)));
    cat_tsp->add_child(std::make_unique<configuration_kv_pair>(
        "Enable covariance transport",
        std::format("{}", m_config.stepping.do_covariance_transport)));

    if (m_config.stepping.do_covariance_transport) {
        auto cat_cov =
            std::make_unique<configuration_category>("Covariance transport");

        cat_cov->add_child(std::make_unique<configuration_kv_pair>(
            "Enable energy loss gradient",
            std::format("{}", m_config.stepping.use_eloss_gradient)));
        cat_cov->add_child(std::make_unique<configuration_kv_pair>(
            "Enable B-field gradient",
            std::format("{}", m_config.stepping.use_field_gradient)));

        cat_tsp->add_child(std::move(cat_cov));
    }

    auto cat_geo = std::make_unique<configuration_category>("Geometry context");

    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::move(cat_nav));
    cat->add_child(std::move(cat_tsp));
    cat->add_child(std::move(cat_geo));

    return cat;
}

}  // namespace traccc::opts
