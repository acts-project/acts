/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"

// Project include(s)
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

namespace traccc::device {

/// Seed track parameter estimation algorithm
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
struct seed_parameter_estimation_algorithm
    : public algorithm<bound_track_parameters_collection_types::buffer(
          const magnetic_field&, const edm::measurement_collection::const_view&,
          const edm::spacepoint_collection::const_view&,
          const edm::seed_collection::const_view&)>,
      public messaging,
      public algorithm_base {

    public:
    /// Constructor for the seed parameter estimation algorithm
    ///
    /// @param config The track parameter estimation configuration
    /// @param mr is the memory resource
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param logger The logger instance to use
    ///
    seed_parameter_estimation_algorithm(
        const track_params_estimation_config& config,
        const traccc::memory_resource& mr, vecmem::copy& copy,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());
    /// Destructor
    ~seed_parameter_estimation_algorithm();

    /// Operator executing the algorithm.
    ///
    /// @param bfield The magnetic field object
    /// @param measurements All measurements of the event
    /// @param spacepoints All spacepoints of the event
    /// @param seeds The reconstructed track seeds of the event
    /// @return A vector of bound track parameters for the seeds
    ///
    output_type operator()(
        const magnetic_field& bfield,
        const edm::measurement_collection::const_view& measurements,
        const edm::spacepoint_collection::const_view& spacepoints,
        const edm::seed_collection::const_view& seeds) const override;

    protected:
    /// @name Function(s) to be implemented by derived classes
    /// @{

    /// Payload for the @c estimate_seed_params_kernel function
    struct estimate_seed_params_kernel_payload {
        /// The number of seeds
        edm::seed_collection::const_view::size_type n_seeds;
        /// The track parameter estimation configuration
        const track_params_estimation_config& config;
        /// The magnetic field object
        const magnetic_field& bfield;
        /// All measurements of the event
        const edm::measurement_collection::const_view& measurements;
        /// All spacepoints of the event
        const edm::spacepoint_collection::const_view& spacepoints;
        /// The reconstructed track seeds of the event
        const edm::seed_collection::const_view& seeds;
        /// The output buffer for the bound track parameters
        bound_track_parameters_collection_types::view& params;
    };

    /// Seed parameter estimation kernel launcher
    ///
    /// @param payload The payload for the kernel
    ///
    virtual void estimate_seed_params_kernel(
        const struct estimate_seed_params_kernel_payload& payload) const = 0;

    /// @}

    private:
    /// Internal data type
    struct data;
    /// Pointer to internal data
    std::unique_ptr<data> m_data;

};  // struct seed_parameter_estimation_algorithm

}  // namespace traccc::device
