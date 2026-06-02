// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/FreeTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"

namespace {
[[deprecated(
    "This header is deprecated. Please include the specific track parameters "
    "header instead.")]]
constexpr static int Acts_EventData_TrackParameters_hpp_is_deprecated = 0;
constexpr static int please_dont_use_Acts_EventData_TrackParameters_hpp =
    Acts_EventData_TrackParameters_hpp_is_deprecated;
}  // namespace
