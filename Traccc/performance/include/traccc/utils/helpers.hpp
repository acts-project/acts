/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <string>
#include <string_view>

namespace traccc::plot_helpers {

/// @brief Nested binning struct for booking plots
struct binning {

    /// Constructor with default arguments
    binning(std::string_view b_title = "", int bins = 0, float b_min = 0.f,
            float b_max = 0.f)
        : title(b_title), n_bins(bins), min(b_min), max(b_max) {}

    std::string title;  ///< title to be displayed
    int n_bins;         ///< number of bins
    float min;          ///< minimum value
    float max;          ///< maximum value
};

}  // namespace traccc::plot_helpers
