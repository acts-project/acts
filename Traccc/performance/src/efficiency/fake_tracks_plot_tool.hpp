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

// This is a duplicate of file "duplication_plot_tool.hpp", adapted to "fake
// tracks".

namespace traccc {

// Tools for creating "fake tracks" plots to show the proportion of fake tracks
// in the reconstruction data.
//
class fake_tracks_plot_tool {
    public:
    /// @brief The nested configuration struct
    struct config {
        std::map<std::string, plot_helpers::binning> var_binning = {
            {"Eta", plot_helpers::binning("#eta", 40, -4.f, 4.f)},
            {"Phi", plot_helpers::binning("#phi", 100, -3.15f, 3.15f)},
            {"Pt", plot_helpers::binning("pT [GeV/c]", 40, 0.f, 100.f)},
            {"Num", plot_helpers::binning("N", 30, -0.5f, 29.5f)}};
    };

    /// @brief Nested Cache struct
    struct fake_tracks_plot_cache {
#ifdef TRACCC_HAVE_ROOT
        std::unique_ptr<TProfile>
            n_fake_vs_pT;  ///< Number of fake tracks vs pT
        std::unique_ptr<TProfile>
            n_fake_vs_eta;  ///< Number of fake tracks vs eta
        std::unique_ptr<TProfile>
            n_fake_vs_phi;  ///< Number of fake tracks vs phi
        std::unique_ptr<TEfficiency>
            fake_rate_vs_pT;  ///< Tracking fake tracks rate vs pT
        std::unique_ptr<TEfficiency>
            fake_rate_vs_eta;  ///< Tracking fake tracks rate vs eta
        std::unique_ptr<TEfficiency>
            fake_rate_vs_phi;  ///< Tracking fake tracks rate vs phi
#endif                         // TRACCC_HAVE_ROOT
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    fake_tracks_plot_tool(const config& cfg) : m_cfg(cfg) {}

    /// @brief book the fake tracks plots
    ///
    /// @param fakePlotCache the cache for fake tracks plots
    void book(std::string_view name, fake_tracks_plot_cache& cache) const {

        plot_helpers::binning b_pt = m_cfg.var_binning.at("Pt");
        plot_helpers::binning b_eta = m_cfg.var_binning.at("Eta");
        plot_helpers::binning b_phi = m_cfg.var_binning.at("Phi");
        plot_helpers::binning b_num = m_cfg.var_binning.at("Num");

        // Avoid unused variable warnings when building the code without ROOT.
        (void)name;
        (void)cache;

#ifdef TRACCC_HAVE_ROOT
        // fake tracks rate vs pT
        cache.fake_rate_vs_pT = plot_helpers::book_eff(
            TString(name) + "_fakeTracksRate_vs_pT",
            "Fake tracks rate;pT [GeV/c];Fake tracks rate", b_pt);
        // fake tracks rate vs eta
        cache.fake_rate_vs_eta = plot_helpers::book_eff(
            TString(name) + "_fakeTracksRate_vs_eta",
            "Fake tracks rate;#eta;Fake tracks rate", b_eta);
        // fake tracks rate vs phi
        cache.fake_rate_vs_phi = plot_helpers::book_eff(
            TString(name) + "_fakeTracksRate_vs_phi",
            "Fake tracks rate;#phi;Fake tracks rate", b_phi);

        // fake tracks number vs pT
        cache.n_fake_vs_pT = plot_helpers::book_prof(
            TString(name) + "_nFakeTracks_vs_pT",
            "Averaged number of fake tracks per particle", b_pt, b_num);
        // fake tracks number vs eta
        cache.n_fake_vs_eta = plot_helpers::book_prof(
            TString(name) + "_nFakeTracks_vs_eta",
            "Averaged number of fake tracks per particle", b_eta, b_num);
        // fake tracks number vs phi
        cache.n_fake_vs_phi = plot_helpers::book_prof(
            TString(name) + "_nFakeTracks_vs_phi",
            "Averaged number of fake tracks per particle", b_phi, b_num);
#endif  // TRACCC_HAVE_ROOT
    }

    /// @brief fill number of fake tracks for a truth particle seed
    ///
    /// @param fakePlotCache cache object for fake tracks plots
    /// @param truthParticle the truth Particle
    /// @param nDuplicatedTracks the number of fake tracks
    void fill(fake_tracks_plot_cache& cache, const particle& truth_particle,
              size_t n_fake_tracks) const {
        const auto t_phi = vector::phi(truth_particle.momentum);
        const auto t_eta = vector::eta(truth_particle.momentum);
        const auto t_pT = vector::perp(
            vector2{truth_particle.momentum[0], truth_particle.momentum[1]});

        // Avoid unused variable warnings when building the code without ROOT.
        (void)t_phi;
        (void)t_eta;
        (void)t_pT;
        (void)cache;
        (void)n_fake_tracks;

#ifdef TRACCC_HAVE_ROOT
        cache.n_fake_vs_pT->Fill(t_pT, static_cast<double>(n_fake_tracks));
        cache.n_fake_vs_eta->Fill(t_eta, static_cast<double>(n_fake_tracks));
        cache.n_fake_vs_phi->Fill(t_phi, static_cast<double>(n_fake_tracks));
#endif  // TRACCC_HAVE_ROOT
    }

    /// @brief write the fake tracks plots to file
    ///
    /// @param fakePlotCache cache object for fake tracks plots
    void write(const fake_tracks_plot_cache& cache) const {

        // Avoid unused variable warnings when building the code without ROOT.
        (void)cache;

#ifdef TRACCC_HAVE_ROOT
        cache.fake_rate_vs_pT->Write();
        cache.fake_rate_vs_eta->Write();
        cache.fake_rate_vs_phi->Write();
        cache.n_fake_vs_pT->Write();
        cache.n_fake_vs_eta->Write();
        cache.n_fake_vs_phi->Write();
#endif  // TRACCC_HAVE_ROOT
    }

    private:
    config m_cfg;  ///< The Config class
};

}  // namespace traccc
