/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Command line options used to configure performance measurements
class performance : public interface {

    public:
    /// @name Options
    /// @{

    /// Whether to run performance checks
    bool run = false;

    /// @}

    /// Constructor
    performance();

    std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct performance

}  // namespace traccc::opts
