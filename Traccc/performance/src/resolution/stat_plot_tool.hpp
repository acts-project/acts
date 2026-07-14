/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "../utils/helpers.hpp"
#include "traccc/resolution/stat_plot_tool_config.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

namespace traccc {

class stat_plot_tool {

    public:
    /// @brief Nested Cache struct
    struct stat_plot_cache {
#ifdef TRACCC_HAVE_ROOT
        // Histogram for the number of DoFs
        std::unique_ptr<TH1> ndf_hist;
        // Histogram for the chi sqaure
        std::unique_ptr<TH1> chi2_hist;
        // Histogram for the chi2/ndf
        std::unique_ptr<TH1> reduced_chi2_hist;
        // Histogram for the pvalue
        std::unique_ptr<TH1> pval_hist;
        // Histogram for the purity
        std::unique_ptr<TH1> purity_hist;
        // Histogram for the completeness
        std::unique_ptr<TH1> completeness_hist;
        // Histogram for chi2 of filtered states
        std::map<unsigned int, std::unique_ptr<TH1>> chi2_filtered_hist;
        // Histogram for chi2 of smoothed states
        std::map<unsigned int, std::unique_ptr<TH1>> chi2_smoothed_hist;
        // Histogram for p-value of filtered states
        std::map<unsigned int, std::unique_ptr<TH1>> pval_filtered_hist;
        // Histogram for p-value of smoothed states
        std::map<unsigned int, std::unique_ptr<TH1>> pval_smoothed_hist;
#endif  // TRACCC_HAVE_ROOT
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    stat_plot_tool(const stat_plot_tool_config& cfg);

    /// @brief book the statistics plots
    ///
    /// @param cache the cache for statistics plots
    void book(stat_plot_cache& cache) const;

    /// @brief fill the cache
    ///
    /// @tparam track_backend_t the backend type for @c fit_res
    ///
    /// @param cache the cache for statistics plots
    /// @param fit_res fitting information that contains statistics
    template <typename track_backend_t>
    void fill(stat_plot_cache& cache,
              const edm::track<track_backend_t>& fit_res) const;

    /// @brief fill the cache
    ///
    /// @tparam track_state_backend_t the backend type for @c trk_state
    ///
    /// @param cache the cache for statistics plots
    /// @param trk_state track state at local measurements
    template <typename track_state_backend_t>
    void fill(stat_plot_cache& cache,
              const edm::track_state<track_state_backend_t>& trk_state,
              const edm::measurement_collection::host& measurements) const;

    /// @brief fill the cache
    /// @param cache the cache for statistics plots
    /// @param purity the track purity
    /// @param completeness the track completeness
    void fill(stat_plot_cache& cache, double purity, double completeness) const;

    /// @brief write the statistics plots into ROOT
    ///
    /// @param cache the cache for statistics plots
    void write(const stat_plot_cache& cache) const;

    private:
    stat_plot_tool_config m_cfg;  ///< The Config class
};

}  // namespace traccc

// Include the implementation.
#include "stat_plot_tool.ipp"
