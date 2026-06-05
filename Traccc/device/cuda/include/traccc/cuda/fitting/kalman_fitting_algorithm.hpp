/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/cuda/utils/stream.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/fitting/fitting_config.hpp"
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

/// Kalman filter based track fitting algorithm
class kalman_fitting_algorithm
    : public algorithm<edm::track_container<default_algebra>::buffer(
          const detector_buffer&, const magnetic_field&,
          const edm::track_container<default_algebra>::const_view&)>,
      public messaging {

    public:
    /// Configuration type
    using config_type = fitting_config;

    /// Constructor with the algorithm's configuration
    ///
    /// @param config The configuration object
    /// @param mr     The memory resource(s) used by the algorithm
    /// @param copy   The copy object used by the algorithm
    /// @param str    The CUDA stream used by the algorithm
    /// @param logger The logger used by the algorithm
    ///
    kalman_fitting_algorithm(
        const config_type& config, const traccc::memory_resource& mr,
        vecmem::copy& copy, stream& str,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Execute the algorithm
    ///
    /// @param det             The detector object
    /// @param field           The magnetic field object
    /// @param track_candidates All track candidates to fit
    ///
    /// @return A container of the fitted track states
    ///
    output_type operator()(
        const detector_buffer& det, const magnetic_field& field,
        const edm::track_container<default_algebra>::const_view&
            track_candidates) const override;

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

};  // class kalman_fitting_algorithm

}  // namespace traccc::cuda
