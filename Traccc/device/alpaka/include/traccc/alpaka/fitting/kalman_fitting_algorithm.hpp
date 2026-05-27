/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/alpaka/utils/queue.hpp"

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

namespace traccc::alpaka {

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
    /// @param q      The Alpaka queue used by the algorithm
    /// @param logger The logger used by the algorithm
    ///
    kalman_fitting_algorithm(
        const config_type& config, const traccc::memory_resource& mr,
        vecmem::copy& copy, queue& q,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Execute the algorithm
    ///
    /// @param det             The detector object
    /// @param bfield          The magnetic field object
    /// @param track_candidates All track candidates to fit
    ///
    /// @return A container of the fitted track states
    ///
    output_type operator()(
        const detector_buffer& det, const magnetic_field& bfield,
        const edm::track_container<default_algebra>::const_view&
            track_candidates) const override;

    private:
    /// Algorithm configuration
    config_type m_config;
    /// Memory resource used by the algorithm
    traccc::memory_resource m_mr;
    /// Copy object used by the algorithm
    std::reference_wrapper<vecmem::copy> m_copy;
    /// Queue wrapper
    std::reference_wrapper<queue> m_queue;

};  // class kalman_fitting_algorithm

}  // namespace traccc::alpaka
