/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/stream.hpp"

// Project include(s).
#include "traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

namespace traccc::cuda {

/// Evicts tracks that seem to be duplicates or fakes. This algorithm takes a
/// greedy approach in the sense that it will remove the track which looks "most
/// duplicate/fake"
class greedy_ambiguity_resolution_algorithm
    : public algorithm<edm::track_container<default_algebra>::buffer(
          const edm::track_container<default_algebra>::const_view&)>,
      public messaging {

    public:
    using config_type = ambiguity_resolution_config;

    /// Constructor for the greedy ambiguity resolution algorithm
    ///
    /// @param cfg The configuration
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param str The CUDA stream to perform the operations in
    ///
    greedy_ambiguity_resolution_algorithm(
        const config_type& cfg, const traccc::memory_resource& mr,
        vecmem::copy& copy, stream& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Run the algorithm
    ///
    /// @param tracks the container view of found patterns
    /// @return the container without ambiguous tracks
    output_type operator()(
        const edm::track_container<default_algebra>::const_view& tracks)
        const override;

    /// Get configuration
    config_type& get_config() { return m_config; }

    private:
    /// Algorithm configuration
    config_type m_config;
    /// The memory resource to use
    traccc::memory_resource m_mr;
    /// The copy object to use
    std::reference_wrapper<vecmem::copy> m_copy;
    /// The CUDA stream to use
    std::reference_wrapper<stream> m_stream;
    /// Warp size of the GPU being used
    unsigned int m_warp_size;
};

}  // namespace traccc::cuda
