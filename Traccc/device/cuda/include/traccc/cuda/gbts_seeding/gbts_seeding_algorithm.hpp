/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/cuda/utils/stream.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// GBTS include(s)
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <memory>

namespace traccc::cuda {

/// Main algorithm for performing GBTS on an NVIDIA GPU
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class gbts_seeding_algorithm
    : public algorithm<edm::seed_collection::buffer(
          const edm::spacepoint_collection::const_view&,
          const edm::measurement_collection::const_view&)>,
      public messaging {

    public:
    /// Constructor for the seed finding algorithm
    ///
    /// @param str The CUDA stream to perform the operations in
    ///
    gbts_seeding_algorithm(
        const gbts_seedfinder_config& cfg, const traccc::memory_resource& mr,
        vecmem::copy& copy, stream& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Operator executing the algorithm.
    ///
    /// @param spacepoints is a view of all spacepoints in the event
    /// @param measurements is a view of all measurements in the event
    /// @return the buffer of track seeds reconstructed from the spacepoints
    ///
    output_type operator()(
        const edm::spacepoint_collection::const_view& spacepoints,
        const edm::measurement_collection::const_view& measurements) const;

    private:
    gbts_seedfinder_config m_config;
    /// The memory resource(s) to use
    traccc::memory_resource m_mr;
    /// The copy object to use
    std::reference_wrapper<vecmem::copy> m_copy;
    /// The CUDA stream to use
    std::reference_wrapper<stream> m_stream;
};

}  // namespace traccc::cuda
