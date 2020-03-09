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
  /// @tparam source_link_t The type of source link
  ///
  /// @param trajectory The multitrajectory object
  /// @param entryIndex The entry index for a track
  ///
  /// @return The resulting
  template <typename source_link_t>
  bool operator()(const MultiTrajectory<source_link_t>& /*trajectory*/,
                  size_t /*entryIndex*/) const {
    return false;
  }
};

}  // namespace Acts
