/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// This is pretty stupid, but "math.h" from either thrust or Vc
// (which are both included as "system includes") confuses GCC
// in the build setup that we get for this source file. Which needs
// to be worked around by including <cmath> early, before one of
// those files (don't know which...) would've been parsed by GCC.
#include <cmath>

// Local include(s).
#include "triplet_fitting_test.hpp"

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#endif  // TRACCC_HAVE_ROOT

// System include(s).
#include <iostream>
#include <memory>
#include <stdexcept>

namespace traccc {

void TripletFittingTests::pull_value_tests(
    std::string_view file_name,
    const std::vector<std::string>& hist_names) const {

    // Avoid unused variable warnings when building the code without ROOT.
    (void)file_name;
    (void)hist_names;

#ifdef TRACCC_HAVE_ROOT
    // Open the file with the histograms.
    std::unique_ptr<TFile> ifile(TFile::Open(file_name.data(), "READ"));
    if ((!ifile) || ifile->IsZombie()) {
        throw std::runtime_error(std::string("Could not open file \"") +
                                 file_name.data() + "\"");
    }

    // Process every (requested) histogram from the file.
    for (const std::string& hname : hist_names) {

        // Access the histogram.
        TH1* pull_dist = dynamic_cast<TH1*>(ifile->Get(hname.c_str()));
        if (!pull_dist) {
            throw std::runtime_error("Could not access histogram \"" + hname +
                                     "\" in file \"" + file_name.data() + "\"");
        }

        // Function used for the fit.
        TF1 gaus{"gaus", "gaus", -5, 5};
        double fit_par[3];

        // Set the mean seed to 0
        gaus.SetParameters(1, 0.);
        gaus.SetParLimits(1, -5., 5.);  // -1., 1.
        // Set the standard deviation seed to 1
        gaus.SetParameters(2, 1.0);
        gaus.SetParLimits(2, 0.1, 5.);  // 0.5, 2.

        auto res = pull_dist->Fit("gaus", "Q0S");

        gaus.GetParameters(&fit_par[0]);

        // Mean check
        EXPECT_NEAR(fit_par[1], 0, 0.05) << hname << " mean value error";

        // Sigma check
        EXPECT_NEAR(fit_par[2], 1, 0.1) << hname << " sigma value error";
    }
#else
    std::cout << "Pull value tests not performed without ROOT" << std::endl;
#endif  // TRACCC_HAVE_ROOT
}

void TripletFittingTests::p_value_tests(
    [[maybe_unused]] std::string_view file_name) const {

#ifdef TRACCC_HAVE_ROOT
    // Open the file with the histograms.
    std::unique_ptr<TFile> ifile(TFile::Open(file_name.data(), "READ"));
    if ((!ifile) || ifile->IsZombie()) {
        throw std::runtime_error(std::string("Could not open file \"") +
                                 file_name.data() + "\"");
    }

    // Access the pvalue histogram.
    TH1* pval_hist = dynamic_cast<TH1*>(ifile->Get("pval"));
    if (!pval_hist) {

        throw std::runtime_error("Could not access p-value histogram in file " +
                                 std::string(file_name.data()) + "\"");
    }

    // Contant function used for the fit.
    TF1 constant_func{"constant", "[0]", 0.f, 1.f};

    // Set the seed to the number of data points divided by the number of bins
    double exp_par = static_cast<double>(pval_hist->GetEntries()) /
                     static_cast<double>(pval_hist->GetNbinsX());

    constant_func.SetParameters(0, exp_par);

    TFitResultPtr res = pval_hist->Fit("constant", "Q0S");

    auto fit_par = constant_func.GetParameters();

    // Some sanity checks
    EXPECT_TRUE(fit_par[0] >= 1.f);
    EXPECT_TRUE(exp_par >= 1.f);

    EXPECT_NEAR(pval_hist->GetMean(), 0.5f, 0.05f);

    // Fitted constant check (5% errror)
    EXPECT_NEAR(fit_par[0], exp_par, exp_par * 0.05f)
        << " fitted constant does not make sense";

    // Make sure the pvalue from constant fit is higher than 0.01
    EXPECT_GE(res->Prob(), 0.01f) << " Poor constant fitting quality";

#else
    std::cout << "P-value tests not performed without ROOT" << std::endl;
#endif  // TRACCC_HAVE_ROOT
}

void TripletFittingTests::ndf_tests(
    const edm::track_collection<default_algebra>::host::const_proxy_type& track,
    const edm::track_state_collection<default_algebra>::host& track_states,
    const edm::measurement_collection::host& measurements) {

    scalar dim_sum = 0;
    std::size_t n_effective_states = 0;

    for (const auto& [type, index] : track.constituent_links()) {

        if (type == edm::track_constituent_link::track_state) {

            auto state = track_states.at(index);

            // Getting rid of the requirement for state
            // to be smoothed, since it does not apply
            // for the Triplet fit
            if (!state.is_hole()) {

                dim_sum += static_cast<scalar>(
                    measurements.at(state.measurement_index()).dimensions());
                n_effective_states++;
            }

        } else if (type == edm::track_constituent_link::measurement) {
            dim_sum += static_cast<scalar>(measurements.at(index).dimensions());
        } else {
            GTEST_FAIL();
        }
    }

    // Check if the number of degree of freedoms is equal to (the sum of
    // measurement dimensions - 5)
    ASSERT_FLOAT_EQ(static_cast<float>(track.ndf()),
                    static_cast<float>(dim_sum) - 5.f);

    // The number of track states is supposed to be eqaul to the number
    // of measurements unless KF failes in the middle of propagation
    if (n_effective_states == track.constituent_links().size()) {
        n_success++;
    }
}

std::size_t TripletFittingTests::count_successfully_fitted_tracks(
    const edm::track_collection<default_algebra>::host& tracks) const {

    const std::size_t n_tracks = tracks.size();
    std::size_t n_fitted_tracks = 0u;

    for (std::size_t i = 0; i < n_tracks; ++i) {
        if (tracks.at(i).fit_outcome() == track_fit_outcome::SUCCESS) {
            n_fitted_tracks++;
        }
    }

    return n_fitted_tracks;
}

}  // namespace traccc
