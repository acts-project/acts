/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "res_plot_tool.hpp"

// Detray include(s).
#include <detray/test/utils/statistics.hpp>

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <TF1.h>
#include <TGraphErrors.h>
#endif  // TRACCC_HAVE_ROOT

namespace traccc {

res_plot_tool::res_plot_tool(const res_plot_tool_config& cfg) : m_cfg(cfg) {}

void res_plot_tool::book(res_plot_cache& cache) const {

    plot_helpers::binning b_pull = m_cfg.var_binning.at("pull");
    plot_helpers::binning b_eta = m_cfg.var_binning.at("Eta");
    plot_helpers::binning b_pT = m_cfg.var_binning.at("Pt");

    // Avoid unused variable warnings when building the code without ROOT.
    (void)cache;
    (void)b_pull;
    (void)b_eta;
    (void)b_pT;

    for (std::size_t idx = 0; idx < m_cfg.param_names.size(); idx++) {
        std::string par_name = m_cfg.param_names.at(idx);

        // Binning for residual is parameter dependent
        std::string par_residual = "residual_" + par_name;
        plot_helpers::binning b_residual = m_cfg.var_binning.at(par_residual);

        // Avoid unused variable warnings when building the code without ROOT.
        (void)b_residual;

#ifdef TRACCC_HAVE_ROOT
        // residual distributions
        cache.residuals[par_name] = plot_helpers::book_histo(
            Form("res_%s", par_name.c_str()),
            Form("Residual of %s", par_name.c_str()), b_residual);

        // pull distritutions
        cache.pulls[par_name] = plot_helpers::book_histo(
            Form("pull_%s", par_name.c_str()),
            Form("Pull of %s", par_name.c_str()), b_pull);

        // resolution vs. eta plots
        cache.resolutions_eta[par_name] = plot_helpers::book_histo(
            Form("resolution_%s_vs_eta", par_name.c_str()),
            Form("Resolution of %s vs. eta", par_name.c_str()), b_eta);

        // resolution vs. pT plots
        cache.resolutions_pT[par_name] = plot_helpers::book_histo(
            Form("resolution_%s_vs_pT", par_name.c_str()),
            Form("Resolution of %s vs. pT", par_name.c_str()), b_pT);

        // residuals vs. eta plots
        cache.residuals_eta[par_name] = plot_helpers::book_histo(
            Form("residual_%s_vs_eta", par_name.c_str()),
            Form("Residual of %s vs. eta", par_name.c_str()), b_eta,
            b_residual);

        // residuals vs. pT plots
        cache.residuals_pT[par_name] = plot_helpers::book_histo(
            Form("residual_%s_vs_pT", par_name.c_str()),
            Form("Residual of %s vs. pT", par_name.c_str()), b_pT, b_residual);

        for (int i = 0; i < b_eta.n_bins; i++) {
            cache.residuals_per_eta[par_name][static_cast<unsigned long>(i)] =
                plot_helpers::book_histo(
                    Form("residual_%s_at_%d_th_eta", par_name.c_str(), i),
                    Form("Residual of %s at %d th eta", par_name.c_str(), i),
                    b_residual);
        }

        for (int i = 0; i < b_pT.n_bins; i++) {
            cache.residuals_per_pT[par_name][static_cast<unsigned long>(i)] =
                plot_helpers::book_histo(
                    Form("residual_%s_at_%d_th_pT", par_name.c_str(), i),
                    Form("Residual of %s at %d th pT", par_name.c_str(), i),
                    b_residual);
        }
#endif  // TRACCC_HAVE_ROOT
    }
}

void res_plot_tool::fill(res_plot_cache& cache,
                         const bound_track_parameters<>& truth_param,
                         const bound_track_parameters<>& fit_param,
                         const particle& ptc) const {

    // Find index of eta and pT for resolution histogram
    const scalar eta = vector::eta(ptc.momentum);
    const scalar pT = std::hypot(ptc.momentum[0], ptc.momentum[1]);

    // Avoid unused variable warnings when building the code without ROOT.
    (void)cache;
    (void)eta;
    (void)pT;

    for (std::size_t idx = 0; idx < m_cfg.param_names.size(); idx++) {
        std::string par_name = m_cfg.param_names.at(idx);

        scalar residual = 0.f;
        scalar pull = 0.f;

        if (idx < e_bound_size) {
            residual = fit_param[idx] - truth_param[idx];
            pull = residual /
                   std::sqrt(getter::element(fit_param.covariance(), idx, idx));
        } else if (par_name == "qopT") {
            residual = fit_param.qopT() - truth_param.qopT();
        } else if (par_name == "qopz") {
            residual = fit_param.qopz() - truth_param.qopz();
        }

        // Avoid unused variable warnings when building the code without ROOT.
        (void)residual;
        (void)pull;

#ifdef TRACCC_HAVE_ROOT
        const auto eta_idx =
            std::clamp(cache.resolutions_eta[par_name]->FindBin(eta) - 1, 0,
                       cache.resolutions_eta[par_name]->GetNbinsX() - 1);
        const auto pT_idx =
            std::clamp(cache.resolutions_pT[par_name]->FindBin(pT) - 1, 0,
                       cache.resolutions_pT[par_name]->GetNbinsX() - 1);

        cache.residuals.at(par_name)->Fill(residual);
        if (idx < e_bound_size) {
            cache.pulls.at(par_name)->Fill(pull);
        }
        cache.residuals_eta.at(par_name)->Fill(eta, residual);
        cache.residuals_pT.at(par_name)->Fill(pT, residual);
        cache.residuals_per_eta.at(par_name)
            .at(static_cast<std::size_t>(eta_idx))
            ->Fill(residual);
        cache.residuals_per_pT.at(par_name)
            .at(static_cast<std::size_t>(pT_idx))
            ->Fill(residual);
#endif  // TRACCC_HAVE_ROOT
    }
}

void res_plot_tool::write(res_plot_cache& cache) const {

    // Avoid unused variable warnings when building the code without ROOT.
    (void)cache;

#ifdef TRACCC_HAVE_ROOT
    for (const auto& residual : cache.residuals) {
        residual.second->Write();
    }
    for (const auto& pull : cache.pulls) {
        pull.second->Write();
    }
    for (const auto& residual : cache.residuals_eta) {
        residual.second->Write();
    }
    for (const auto& residual : cache.residuals_pT) {
        residual.second->Write();
    }

    auto create_resolution_graph = [](std::shared_ptr<TH1> H, auto residuals,
                                      const char* x_axis_title,
                                      const char* y_axis_title) {
        const auto n_bins = H->GetNbinsX();

        std::vector<float> sigmas;

        for (int i = 0; i < n_bins; i++) {
            auto& data = residuals[static_cast<std::size_t>(i)];

            // When there is no entry
            if (data->GetEntries() == 0) {
                H->SetBinContent(i, 0);
                sigmas.push_back(0);
                continue;
            }

            // Function used for the fit.
            TF1 gaus{"gaus", "gaus", -1., 1.};
            double fit_par[3];
            auto res = data->Fit("gaus", "Q0S");
            gaus.GetParameters(&fit_par[0]);
            H->SetBinContent(i + 1, fit_par[2]);
            sigmas.push_back(static_cast<float>(gaus.GetParError(2)));
        }

        std::unique_ptr<TGraphErrors> G =
            std::make_unique<TGraphErrors>(H.get());

        for (int i = 0; i < n_bins; i++) {
            G->SetPointError(i, H->GetBinWidth(i) / 2.f,
                             sigmas[static_cast<std::size_t>(i)]);
        }
        G->SetName(H->GetName());
        G->SetMinimum(1e-5);
        G->SetMaximum(1e-1);
        G->GetHistogram()->GetYaxis()->SetMaxDigits(2);
        G->GetHistogram()->GetXaxis()->SetTitle(x_axis_title);
        G->GetHistogram()->GetYaxis()->SetTitle(y_axis_title);

        return G;
    };

    for (const auto& resolution : cache.resolutions_eta) {
        const auto& par_name = resolution.first;
        auto& H = resolution.second;

        plot_helpers::binning binning =
            m_cfg.var_binning.at("resolution_" + par_name);

        auto G =
            create_resolution_graph(H, cache.residuals_per_eta.at(par_name),
                                    "#eta", binning.title.c_str());

        G->Write();
    }

    for (const auto& resolution : cache.resolutions_pT) {
        const auto& par_name = resolution.first;
        auto& H = resolution.second;

        plot_helpers::binning binning =
            m_cfg.var_binning.at("resolution_" + par_name);

        auto G =
            create_resolution_graph(H, cache.residuals_per_pT.at(par_name),
                                    "p_{T} [GeV/c]", binning.title.c_str());

        G->Write();
    }

#endif  // TRACCC_HAVE_ROOT
}

}  // namespace traccc
