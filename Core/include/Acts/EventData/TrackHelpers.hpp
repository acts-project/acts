// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// This header is deprecated and will be removed in the future.
// Please use the new header instead:
#include "Acts/Utilities/TrackHelpers.hpp"

namespace Acts {

// Funny header deprecation strategy
namespace {
[[deprecated(
    "This header is deprecated - use "
    "Acts/Utilities/TrackHelpers.hpp")]] constexpr static int
    utilities_trackhelpers_hpp_is_deprecated = 0;
constexpr static int please_dont_use_utilities_trackhelpers_hpp =
    utilities_trackhelpers_hpp_is_deprecated;
}  // namespace

}  // namespace Acts
