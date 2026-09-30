// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
