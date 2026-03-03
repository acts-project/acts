// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include "Mille/MilleRecord.h"

namespace ActsPlugins::ActsToMille {

/// The MilleRecord is the primary interface for
/// writing out alignment fit inputs
using Mille::MilleRecord;

/// Placeholder method to test header
/// discovery and linkage.
/// Using TrackAlignmentState as a potential
/// candidate for a future (internal) interface.
void dumpToMille(const ActsAlignment::detail::TrackAlignmentState&,
                 MilleRecord* record);

}  // namespace ActsPlugins::ActsToMille
