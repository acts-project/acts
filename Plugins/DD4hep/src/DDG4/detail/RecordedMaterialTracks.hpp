// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialInteraction.hpp"

#include <map>

namespace ActsPlugins::detail {

/// @brief Material tracks recorded in one event
///
/// Attached to the DDG4 event as an extension by `ActsMaterialRecorder`, and
/// written at the end of the event by `ActsMaterialTrackWriter`.
struct RecordedMaterialTracks {
  /// The material tracks, ordered by Geant4 track ID
  std::map<int, Acts::RecordedMaterialTrack> tracks;
};

}  // namespace ActsPlugins::detail
