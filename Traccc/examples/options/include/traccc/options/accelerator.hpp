/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts {

/// Option(s) for accelerator usage
class accelerator : public interface {

    public:
    /// @name Options
    /// @{

    /// Whether to compare the accelerator code's output with that of the CPU
    bool compare_with_cpu = false;
    /// Whether GPU texture memory should be used
    bool use_gpu_texture_memory = false;

    /// @}

    /// Constructor
    accelerator();

    std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct accelerator

}  // namespace traccc::opts
