/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Options for the cell clusterization algorithm(s)
class clusterization
    : public interface,
      public config_provider<clustering_config>,
      public config_provider<host::clusterization_algorithm::config_type> {

    public:
    /// Constructor
    clusterization();

    /// Configuration conversion
    operator clustering_config() const override;
    operator host::clusterization_algorithm::config_type() const override;

    std::unique_ptr<configuration_printable> as_printable() const override;

    private:
    /// Internal configuration object
    clustering_config m_config;
};  // class clusterization

}  // namespace traccc::opts
