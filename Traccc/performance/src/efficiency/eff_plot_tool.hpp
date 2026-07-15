/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../utils/helpers.hpp"

// Project include(s).
#include "traccc/edm/particle.hpp"

// System include(s).
#include <string_view>

namespace traccc {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class eff_plot_tool {
    public:
    /// @brief The nested configuration struct
    struct config {
        std::map<std::string, plot_helpers::binning> var_binning = {
            {"Eta", plot_helpers::binning("#eta", 40, -4.f, 4.f)},
            {"Phi", plot_helpers::binning("#phi", 100, -3.15f, 3.15f)},
            {"Pt", plot_helpers::binning("pT [GeV/c]", 40, 0.f, 100.f)}};
    };

    /// @brief Nested Cache struct
    struct eff_plot_cache {
#ifdef TRACCC_HAVE_ROOT
        std::unique_ptr<TEfficiency>
            track_eff_vs_pT;  ///< Tracking efficiency vs pT
        std::unique_ptr<TEfficiency>
            track_eff_vs_eta;  ///< Tracking efficiency vs eta
        std::unique_ptr<TEfficiency>
            track_eff_vs_phi;  ///< Tracking efficiency vs phi
#endif                         // TRACCC_HAVE_ROOT
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param lvl Message level declaration
    eff_plot_tool(const config& cfg) : m_cfg(cfg) {}

    /// @brief book the efficiency plots
    ///
    /// @param effPlotCache the cache for efficiency plots
    void book(std::string_view name, eff_plot_cache& cache) const {

        plot_helpers::binning b_phi = m_cfg.var_binning.at("Phi");
        plot_helpers::binning b_eta = m_cfg.var_binning.at("Eta");
        plot_helpers::binning b_pt = m_cfg.var_binning.at("Pt");

        // Avoid unused variable warnings when building the code without ROOT.
        (void)name;
        (void)cache;

#ifdef TRACCC_HAVE_ROOT
        std::string header;

        if (name == "finding") {
            header = "Track finding efficiency";
        } else if (name == "seeding") {
            header = "Seed finding efficiency";
        } else {
            header = "Tracking efficiency";
        }

        // efficiency vs pT
        cache.track_eff_vs_pT = plot_helpers::book_eff(
            TString(name) + "_trackeff_vs_pT",
            header + ";Truth pT [GeV/c];Efficiency", b_pt);
        // efficiency vs eta
        cache.track_eff_vs_eta =
            plot_helpers::book_eff(TString(name) + "_trackeff_vs_eta",
                                   header + ";Truth #eta;Efficiency", b_eta);
        // efficiency vs phi
        cache.track_eff_vs_phi =
            plot_helpers::book_eff(TString(name) + "_trackeff_vs_phi",
                                   header + ";Truth #phi;Efficiency", b_phi);
#endif  // TRACCC_HAVE_ROOT
    }

    /// @brief fill efficiency plots
    ///
    /// @param effPlotCache cache object for efficiency plots
    /// @param truthParticle the truth Particle
    /// @param status the reconstruction status
    void fill(eff_plot_cache& cache, const particle& truth_particle,
              bool status) const {

        const auto t_phi = vector::phi(truth_particle.momentum);
        const auto t_eta = vector::eta(truth_particle.momentum);
        const auto t_pT = vector::perp(
            vector2{truth_particle.momentum[0], truth_particle.momentum[1]});

        // Avoid unused variable warnings when building the code without ROOT.
        (void)t_phi;
        (void)t_eta;
        (void)t_pT;
        (void)cache;
        (void)status;

#ifdef TRACCC_HAVE_ROOT
        cache.track_eff_vs_pT->Fill(status, t_pT);
        cache.track_eff_vs_eta->Fill(status, t_eta);
        cache.track_eff_vs_phi->Fill(status, t_phi);
#endif  // TRACCC_HAVE_ROOT
    }

    /// @brief write the efficiency plots to file
    ///
    /// @param effPlotCache cache object for efficiency plots
    void write(const eff_plot_cache& cache) const {

        // Avoid unused variable warnings when building the code without ROOT.
        (void)cache;

#ifdef TRACCC_HAVE_ROOT
        cache.track_eff_vs_pT->Write();
        cache.track_eff_vs_eta->Write();
        cache.track_eff_vs_phi->Write();
#endif  // TRACCC_HAVE_ROOT
    }

    private:
    config m_cfg;  ///< The Config class
};

}  // namespace traccc
