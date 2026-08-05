/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_container.hpp"
#include "traccc/utils/logging.hpp"

namespace traccc::details {

/// Print statistics for the fitted tracks.
///
/// @param tracks The fitted tracks to print statistics for
/// @param log    The logger to use for outputting the statistics
///
void print_fitted_tracks_statistics(
    const edm::track_container<default_algebra>::host& tracks,
    const Logger& log);

}  // namespace traccc::details
