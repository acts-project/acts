/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/seeding/device/seed_parameter_estimation_algorithm.hpp"

namespace traccc::alpaka {

/// Algorithm for estimating track parameters for seeds using Alpaka
struct seed_parameter_estimation_algorithm
    : public device::seed_parameter_estimation_algorithm,
      public alpaka::algorithm_base {

    public:
    /// Constructor for track_params_estimation
    ///
    /// @param config The track parameter estimation configuration
    /// @param mr is a struct of memory resources (shared or host & device)
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param queue The Alpaka queue to work with
    /// @param logger The logger instance to use
    ///
    seed_parameter_estimation_algorithm(
        const track_params_estimation_config& config,
        const traccc::memory_resource& mr, vecmem::copy& copy,
        alpaka::queue& queue,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    private:
    /// @name Function(s) inherited from
    /// @c traccc::device::seed_parameter_estimation_algorithm
    /// @{

    /// Seed parameter estimation kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    void estimate_seed_params_kernel(
        const struct estimate_seed_params_kernel_payload& payload)
        const override;

    /// @}

};  // struct seed_parameter_estimation_algorithm

}  // namespace traccc::alpaka
