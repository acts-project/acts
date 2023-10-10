// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"

#include <vector>

namespace ActsExamples {

/// (Reconstructed) track parameters e.g. close to the vertex.
using TrackParameters = ::Acts::BoundTrackParameters;
/// Container of reconstructed track states for multiple tracks.
using TrackParametersContainer = std::vector<TrackParameters>;

using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer,
                         Acts::VectorMultiTrajectory, std::shared_ptr>;

using ConstTrackContainer =
    Acts::TrackContainer<Acts::ConstVectorTrackContainer,
                         Acts::ConstVectorMultiTrajectory, std::shared_ptr>;

}  // namespace ActsExamples
