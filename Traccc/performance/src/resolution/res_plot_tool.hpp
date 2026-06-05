/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "../utils/helpers.hpp"
#include "traccc/resolution/res_plot_tool_config.hpp"

// Project include(s).
#include "traccc/edm/particle.hpp"
#include "traccc/edm/track_parameters.hpp"

// System include(s).
#include <map>
#include <memory>
#include <string>

namespace traccc {

class res_plot_tool {

    public:
    /// @brief Nested Cache struct
    struct res_plot_cache {
#ifdef TRACCC_HAVE_ROOT
        // Residuals and pulls for parameters
        std::map<std::string, std::unique_ptr<TH1>> residuals;
        std::map<std::string, std::unique_ptr<TH1>> pulls;
        std::map<std::string, std::shared_ptr<TH1>> resolutions_eta;
        std::map<std::string, std::shared_ptr<TH1>> resolutions_pT;
        std::map<std::string, std::unique_ptr<TH2>> residuals_eta;
        std::map<std::string, std::unique_ptr<TH2>> residuals_pT;
        std::map<std::string, std::map<std::size_t, std::shared_ptr<TH1>>>
            residuals_per_eta;
        std::map<std::string, std::map<std::size_t, std::shared_ptr<TH1>>>
            residuals_per_pT;
#endif  // TRACCC_HAVE_ROOT
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    res_plot_tool(const res_plot_tool_config& cfg);

    /// @brief book the resolution plots
    ///
    /// @param cache the cache for resolution plots
    void book(res_plot_cache& cache) const;

    /// @brief fill the cache
    ///
    /// @param cache the cache for resolution plots
    /// @param truth_param truth track parameter
    /// @param fit_param fitted track parameter
    void fill(res_plot_cache& cache,
              const bound_track_parameters<>& truth_param,
              const bound_track_parameters<>& fit_param,
              const particle& ptc) const;

    /// @brief write the resolution plots into ROOT
    ///
    /// @param cache the cache for resolution plots
    void write(res_plot_cache& cache) const;

    private:
    res_plot_tool_config m_cfg;  ///< The Config class
};

}  // namespace traccc
