// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <tuple>

/// @note this is foreseen for the 'Geometry' module

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;
class Detector;

/// Base class for all link implementations that need class structure
class IDelegateImpl {
 public:
  virtual ~IDelegateImpl() {}
};

/// Memory managed delegate to guarantee the lifetime
/// of eventual unterlying delegate memory and the
/// delegate function
///
template <typename deletage_type>
struct ManagedDelegate {
 public:
  deletage_type delegate;
  std::shared_ptr<IDelegateImpl> implementation = nullptr;
};

/// Declare a navigation state updator
///
/// This delegate dispatches the local navigation action
/// to a dedicated struct or function that is optimised for
/// the given environment.
///
/// @param gctx is the current geometry context
/// @param nState is the navigation state to be updated
///
/// @note it relies on the detector volume to be set to the state
using NavigationStateUpdator =
    Delegate<void(const GeometryContext& gctx, NavigationState& nState)>;

/// Memory  managed navigation state updator
using ManagedNavigationStateUpdator = ManagedDelegate<NavigationStateUpdator>;

/// Declare a Detctor Volume Switching delegate
///
/// @param gctx is the current geometry context
/// @param position is the position at the query
/// @param direction is the direction at the query
///
/// @return the new DetectorVolume into which one changes at this switch
using DetectorVolumeLink = Delegate<const DetectorVolume*(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction)>;

/// Memory managed detector volume link
using ManagedDetectorVolumeLink = ManagedDelegate<DetectorVolumeLink>;

/// @brief  Definition of a volume finder, this can be set and optimised at construction
///
/// @param gctx the geometry context of this call
/// @param detector is the detector in which the volume should be found
/// @param position the search position for this associated volume
///
using DetectorVolumeFinder = Delegate<const DetectorVolume*(
    const GeometryContext& gctx, const Detector& detector,
    const Vector3& position)>;
using ManagedDetectorVolumeFinder = ManagedDelegate<DetectorVolumeFinder>;

}  // namespace Experimental
}  // namespace Acts
