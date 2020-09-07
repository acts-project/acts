// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"

namespace ActsExamples {

/// @brief Helper struct that stores a SimVertex object
/// together with std::vector<Acts::BoundTrackParameters>
struct VertexAndTracks {
  // The vertex
  SimVertex vertex;
  // The tracks
  std::vector<Acts::BoundTrackParameters> tracks;
};

}  // namespace ActsExamples
