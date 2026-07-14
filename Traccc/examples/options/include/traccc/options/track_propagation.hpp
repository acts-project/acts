/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/options/details/value_array.hpp"

// Detray include(s).
#include <detray/propagator/propagation_config.hpp>

namespace traccc::opts {

/// Command line options used in the propagation tests
class track_propagation : public interface,
                          public config_provider<detray::propagation::config> {

    public:
    /// Constructor
    track_propagation();

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    /// Configuration provider
    operator detray::propagation::config() const override;

    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    /// @name Options
    /// @{
    /// Propagation configuration object
    detray::propagation::config m_config;
    /// @}

    /// Search window (helper variable)
    value_array<unsigned int, 2> m_search_window = {0u, 0u};

};  // class track_propagation

}  // namespace traccc::opts
