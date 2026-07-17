/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Configuration for track fitting
class track_fitting : public interface, public config_provider<fitting_config> {

    public:
    /// Constructor
    track_fitting();

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    /// Configuration conversion operators
    operator fitting_config() const override;

    /// Return a printable set of options
    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    /// The internal configuration
    fitting_config m_config;
    /// Particle hypothesis (PDG)
    int m_pdg{m_config.ptc_hypothesis.pdg_num()};

};  // class track_fitting

}  // namespace traccc::opts
