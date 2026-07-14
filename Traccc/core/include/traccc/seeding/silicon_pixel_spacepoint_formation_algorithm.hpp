/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>

namespace traccc::host {

/// Algorithm forming space points out of measurements
///
/// This algorithm performs the local-to-global transformation of the 2D
/// measurements made on every detector module, into 3D spacepoint coordinates.
///
class silicon_pixel_spacepoint_formation_algorithm
    : public algorithm<edm::spacepoint_collection::host(
          const host_detector&,
          const edm::measurement_collection::const_view&)>,
      public messaging {

    public:
    /// Output type
    using output_type = edm::spacepoint_collection::host;

    /// Constructor for spacepoint_formation
    ///
    /// @param mr is the memory resource
    ///
    silicon_pixel_spacepoint_formation_algorithm(
        vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Construct spacepoints from 2D silicon pixel measurements
    ///
    /// @param det Detector object
    /// @param measurements A collection of measurements
    /// @return A spacepoint container, with one spacepoint for every
    ///         silicon pixel measurement
    ///
    output_type operator()(const host_detector& det,
                           const edm::measurement_collection::const_view&
                               measurements) const override;

    private:
    /// Memory resource to use for the output container
    std::reference_wrapper<vecmem::memory_resource> m_mr;
};  // class silicon_pixel_spacepoint_formation_algorithm

}  // namespace traccc::host
