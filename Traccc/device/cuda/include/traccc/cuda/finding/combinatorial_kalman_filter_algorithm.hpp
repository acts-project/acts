/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/cuda/utils/stream.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>

namespace traccc::cuda {

/// CKF track finding algorithm
class combinatorial_kalman_filter_algorithm
    : public algorithm<edm::track_container<default_algebra>::buffer(
          const detector_buffer&, const magnetic_field&,
          const edm::measurement_collection::const_view&,
          const bound_track_parameters_collection_types::const_view&)>,
      public messaging {

    public:
    /// Configuration type
    using config_type = finding_config;

    /// Constructor with the algorithm's configuration
    combinatorial_kalman_filter_algorithm(
        const config_type& config, const traccc::memory_resource& mr,
        vecmem::copy& copy, stream& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Execute the algorithm
    ///
    /// @param det          The (default) detector object
    /// @param bfield       The magnetic field object
    /// @param measurements All measurements in an event
    /// @param seeds        All seeds in an event to start the track finding
    ///                     with
    ///
    /// @return A container of the found track candidates
    ///
    output_type operator()(
        const detector_buffer& det, const magnetic_field& bfield,
        const edm::measurement_collection::const_view& measurements,
        const bound_track_parameters_collection_types::const_view& seeds)
        const override;

    private:
    /// Algorithm configuration
    config_type m_config;
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::reference_wrapper<vecmem::copy> m_copy;
    /// The CUDA stream to use
    std::reference_wrapper<stream> m_stream;
    /// Warp size of the GPU being used
    unsigned int m_warp_size;

};  // class combinatorial_kalman_filter_algorithm

}  // namespace traccc::cuda
