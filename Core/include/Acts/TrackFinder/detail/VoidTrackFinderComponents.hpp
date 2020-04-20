// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <limits>

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
///
/// @brief void branch stopper
struct VoidBranchStopper {
  /// @brief Public call mimicking an track branch stopper
  ///
  /// @tparam track_quality_t The type defining track quality
  ///
  /// @param trackQuality The quality object of a track
  ///
  /// @return The resulting
  template <typename track_quality_t>
  bool operator()(const track_quality_t& /*trackQuality*/) const {
    return false;
  }
};

}  // namespace Acts
