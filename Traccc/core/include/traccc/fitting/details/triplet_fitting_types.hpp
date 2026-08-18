/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/fitting/triplet_filter/triplet_fitter.hpp"

// System include(s).
#include <type_traits>

namespace traccc::details {

/// Triplet fitter type for a specific detector and magnetic field type
///
/// @tparam detector_t The detector type to use
/// @tparam bfield_t   The magnetic field type to use
///
using triplet_fitter_t = traccc::triplet_fitter<
    const typename traccc::default_detector::host,
    typename detray::bfield::const_field_t<
        traccc::default_detector::host::scalar_type>::view_t>;

}  // namespace traccc::details
