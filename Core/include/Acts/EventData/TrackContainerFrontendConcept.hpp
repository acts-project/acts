// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackContainerBackendConcept.hpp"
#include "Acts/EventData/Types.hpp"

#include <concepts>

namespace Acts {

template <typename T>
concept TrackContainerFrontend =
    requires() {
      typename T::ReadOnly;

      requires std::same_as<typename T::IndexType, TrackIndexType>;

      TrackContainerBackend<typename T::TrackContainerBackend>;
      CommonMultiTrajectoryBackend<typename T::TrackStateContainerBackend>;

      typename T::TrackProxy;
      typename T::ConstTrackProxy;
      typename T::TrackStateProxy;
      typename T::ConstTrackStateProxy;
    };

}  // namespace Acts
