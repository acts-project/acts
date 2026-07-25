/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/utils/helpers.hpp"

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#endif  // TRACCC_HAVE_ROOT

// System include(s).
#include <memory>
#include <string_view>

namespace traccc::plot_helpers {

#ifdef TRACCC_HAVE_ROOT

/// @brief book a 1D histogram
/// @param histName the name of histogram
/// @param histTitle the title of histogram
/// @param varBinning the binning info of variable
/// @return histogram pointer
///
std::unique_ptr<TH1> book_histo(std::string_view hist_name,
                                std::string_view hist_title,
                                const binning& var_binning);

/// @brief book a 1D histogram
/// @param hist_name the name of histogram
/// @param hist_title the title of histogram
/// @param var_x_binning the binning info of x variable
/// @param var_y_binning the binning info of y variable
/// @return histogram pointer
///
std::unique_ptr<TH2> book_histo(std::string_view hist_name,
                                std::string_view hist_title,
                                const binning& var_x_binning,
                                const binning& var_y_binning);

/// @brief book a 1D efficiency plot
/// @param effName the name of plot
/// @param effTitle the title of plot
/// @param varBinning the binning info of variable
/// @return TEfficiency pointer
///
std::unique_ptr<TEfficiency> book_eff(std::string_view eff_name,
                                      std::string_view eff_title,
                                      const binning& var_binning);

std::unique_ptr<TProfile> book_prof(std::string_view prof_name,
                                    std::string_view prof_title,
                                    const binning& var_x_binning,
                                    const binning& var_y_binning);

#endif  // TRACCC_HAVE_ROOT

}  // namespace traccc::plot_helpers
