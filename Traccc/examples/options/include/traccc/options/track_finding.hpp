/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/options/details/value_array.hpp"

// System include(s).
#include <limits>

namespace traccc::opts {

/// Configuration for track finding
class track_finding : public interface, public config_provider<finding_config> {

    public:
    /// Constructor
    track_finding();

    /// Read/process the command line options
    ///
    /// @param vm The command line options to interpret/read
    ///
    void read(const boost::program_options::variables_map& vm) override;

    /// Configuration conversion operators
    operator finding_config() const override;

    /// Return a printable set of options
    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    /// The internal configuration
    finding_config m_config;
    /// Additional variables which we cannot simply store in the config
    opts::value_array<unsigned int, 2> m_track_candidates_range{
        m_config.min_track_candidates_per_track,
        m_config.max_track_candidates_per_track};
    /// Helper variable for setting the PDG ID of the particle hypothesis
    int m_pdg_number{m_config.ptc_hypothesis.pdg_num()};

};  // class track_finding

}  // namespace traccc::opts
