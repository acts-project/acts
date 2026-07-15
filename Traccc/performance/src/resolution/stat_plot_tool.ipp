/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/prob.hpp"

namespace traccc {

template <typename track_backend_t>
void stat_plot_tool::fill(stat_plot_cache& cache,
                          const edm::track<track_backend_t>& fit_res) const {

    // Avoid unused variable warnings when building the code without ROOT.
    (void)cache;
    (void)fit_res;

#ifdef TRACCC_HAVE_ROOT
    const scalar ndf = fit_res.ndf();
    const scalar chi2 = fit_res.chi2();
    cache.ndf_hist->Fill(ndf);
    cache.chi2_hist->Fill(chi2);
    cache.reduced_chi2_hist->Fill(chi2 / ndf);
    cache.pval_hist->Fill(fit_res.pval());
#endif  // TRACCC_HAVE_ROOT
}

template <typename track_state_backend_t>
void stat_plot_tool::fill(
    stat_plot_cache& cache,
    const edm::track_state<track_state_backend_t>& trk_state,
    const edm::measurement_collection::host& measurements) const {

    // Avoid unused variable warnings when building the code without ROOT.
    (void)cache;
    (void)trk_state;
    (void)measurements;

#ifdef TRACCC_HAVE_ROOT
    const unsigned int D =
        measurements.at(trk_state.measurement_index()).dimensions();
    const scalar filtered_chi2 = trk_state.filtered_chi2();
    const scalar smoothed_chi2 = trk_state.smoothed_chi2();
    cache.chi2_filtered_hist[D]->Fill(filtered_chi2);
    cache.chi2_smoothed_hist[D]->Fill(smoothed_chi2);
    cache.pval_filtered_hist[D]->Fill(
        prob(filtered_chi2, static_cast<traccc::scalar>(D)));
    cache.pval_smoothed_hist[D]->Fill(
        prob(smoothed_chi2, static_cast<traccc::scalar>(D)));
#endif  // TRACCC_HAVE_ROOT
}

}  // namespace traccc
