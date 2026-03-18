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

/// The MilleRecord is Millepede's interface for
/// writing out alignment fit inputs. It can be instantiated
/// using Mille::spawnMilleRecord(desired_file_name),
/// provided by Mille/MilleFactory.h
using Mille::MilleRecord;

/// @brief Dump a Kalman track encoded as a TrackAlignmentState into
/// a Mille record.
/// @param state: Alignment state to dump.
/// @param record: Mille record to write to - should be valid pointer
/// Note: Not very efficient - we have to "un-fit" the kalman track.
/// Used for R&D, recommending the GBL track model (under development)
//  for production use.
void dumpToMille(const ActsAlignment::detail::TrackAlignmentState& state,
                 MilleRecord* record);

}  // namespace ActsPlugins::ActsToMille
