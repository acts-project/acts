/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// Detray include(s).
#include <covfie/core/field.hpp>

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <functional>
#include <tuple>

namespace traccc::host {

/// Kalman filter based track fitting algorithm
class kalman_fitting_algorithm
    : public algorithm<edm::track_container<default_algebra>::host(
          const host_detector&, const magnetic_field&,
          const edm::track_container<default_algebra>::const_view&)>,
      public messaging {

    public:
    /// Configuration type
    using config_type = fitting_config;

    /// Constructor with the algorithm's configuration
    ///
    /// @param config The configuration object
    ///
    explicit kalman_fitting_algorithm(
        const config_type& config, vecmem::memory_resource& mr,
        vecmem::copy& copy,
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
        const host_detector& det, const magnetic_field& bfield,
        const edm::track_container<default_algebra>::const_view&
            track_candidates) const override;

    private:
    /// Algorithm configuration
    config_type m_config;
    /// Memory resource to use in the algorithm
    std::reference_wrapper<vecmem::memory_resource> m_mr;
    /// The copy object to use
    std::reference_wrapper<vecmem::copy> m_copy;
};  // class kalman_fitting_algorithm

}  // namespace traccc::host
