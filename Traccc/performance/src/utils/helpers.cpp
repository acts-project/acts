/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "helpers.hpp"

namespace traccc::plot_helpers {

#ifdef TRACCC_HAVE_ROOT

std::unique_ptr<TH1> book_histo(std::string_view hist_name,
                                std::string_view hist_title,
                                const binning& var_binning) {
    auto result = std::make_unique<TH1F>(
        hist_name.data(),
        (std::string(hist_title) + ";" + var_binning.title + ";Entries")
            .c_str(),
        var_binning.n_bins, var_binning.min, var_binning.max);
    result->Sumw2();
    return result;
}

std::unique_ptr<TH2> book_histo(std::string_view hist_name,
                                std::string_view hist_title,
                                const binning& var_x_binning,
                                const binning& var_y_binning) {

    auto result = std::make_unique<TH2F>(
        hist_name.data(), hist_title.data(), var_x_binning.n_bins,
        var_x_binning.min, var_x_binning.max, var_y_binning.n_bins,
        var_y_binning.min, var_y_binning.max);

    result->GetXaxis()->SetTitle(var_x_binning.title.c_str());
    result->GetYaxis()->SetTitle(var_y_binning.title.c_str());
    result->Sumw2();
    return result;
}

std::unique_ptr<TEfficiency> book_eff(std::string_view eff_name,
                                      std::string_view eff_title,
                                      const binning& var_binning) {
    return std::make_unique<TEfficiency>(eff_name.data(), eff_title.data(),
                                         var_binning.n_bins, var_binning.min,
                                         var_binning.max);
}

std::unique_ptr<TProfile> book_prof(std::string_view prof_name,
                                    std::string_view prof_title,
                                    const binning& var_x_binning,
                                    const binning& var_y_binning) {

    return std::make_unique<TProfile>(
        prof_name.data(),
        (std::string(prof_title) + ";" + var_x_binning.title + ";" +
         var_y_binning.title)
            .c_str(),
        var_x_binning.n_bins, var_x_binning.min, var_x_binning.max,
        var_y_binning.min, var_y_binning.max);
}

#endif  // TRACCC_HAVE_ROOT

}  // namespace traccc::plot_helpers
