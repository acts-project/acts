/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>
#include <memory>

namespace traccc::host {

/// Algorithm sorting the reconstructed measurements in their container
///
/// The track finding algorithm expects measurements belonging to a single
/// detector module to be consecutive in memory. But certain clusterization /
/// measurement creation algorithms may not produce the measurements in such
/// an ordered state. In such cases this algorithm can be used to sort the
/// measurements "correctly" in place.
///
class measurement_sorting_algorithm
    : public algorithm<edm::measurement_collection::host(
          const edm::measurement_collection::const_view&)>,
      public messaging {

    public:
    /// Constructor
    ///
    /// @param mr The memory resource to use for the algorithm
    /// @param logger The logger to use for the algorithm
    ///
    measurement_sorting_algorithm(
        vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Callable operator performing the sorting on a container
    ///
    /// @param measurements The measurements to sort
    ///
    [[nodiscard]] output_type operator()(
        const edm::measurement_collection::const_view& measurements)
        const override;

    private:
    /// The memory resource to use
    std::reference_wrapper<vecmem::memory_resource> m_mr;

};  // class measurement_sorting_algorithm

}  // namespace traccc::host
