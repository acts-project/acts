// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"

#include "Mille/IMilleReader.h"
#include "Mille/MilleDecoder.h"
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
/// @param record: Mille record to write to.
/// Note: Not very efficient - we have to "un-fit" the kalman track.
/// Used for R&D, recommending the GBL track model (under development)
/// for production use.
void dumpToMille(const ActsAlignment::detail::TrackAlignmentState& state,
                 MilleRecord& record);

/// @brief read one record (= track or (constrained) track pair) from
/// a Mille binary into the equivalent matrices of a TrackAlignmentState.
/// Allows to use Mille to collect tracks across multiple events and
/// align them with the ACTS solver, and to validate the outputs of dumpToMille.
/// @param reader: A Mille Reader, connected to a valid input file.
/// @param targetState: The TrackAlignmentState to populate.
/// @param idxedAlignSurfaces: [optional]: Indexed alignment surfaces from the geometry. If passed,
/// the internal `alignedSurfaces` member of the state will be configured to
/// link back to the correct surfaces.
/// @return a ReadResult enum with 3 possible states to indicate the outcome- ok / end-of-file / read-error.
/// The targetState will only be modified if the result is 'ok'.
Mille::MilleDecoder::ReadResult unpackMilleRecord(
    Mille::IMilleReader& reader,
    ActsAlignment::detail::TrackAlignmentState& targetState,
    const std::unordered_map<const Acts::Surface*, std::size_t>&
        idxedAlignSurfaces);

}  // namespace ActsPlugins::ActsToMille
