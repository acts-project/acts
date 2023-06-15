// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

class Surface;

namespace Experimental {

class ISurfaceCandidatesDelegate {
 public:
  virtual ~ISurfaceCandidatesDelegate() = default;

  virtual void update(const GeometryContext& gctx,
                      const NavigationState& nState,
                      NavigationState::SurfaceCandidates& candidates) const = 0;
};

class IDetectorVolumeFinder {
 public:
  virtual ~IDetectorVolumeFinder() = default;

  virtual const DetectorVolume* find(const GeometryContext& gctx,
                                     const NavigationState& nState) const = 0;
};

/// Declare an updator for the local navigation, i.e. the
/// navigation inside a detector volume. This can be called
/// either directly after a volume switch or in order to update
/// within a volume after some progression
///
/// This delegate dispatches the local navigation action
/// to a dedicated struct or function that is optimised for
/// the given environment.
///
/// @param gctx is the current geometry context
/// @param nState [in,out] is the navigation state to be updated
///
/// @note it relies on the detector volume to be set to the state
/// Memory  managed navigation state updator
using SurfaceCandidatesDelegate =
    OwningDelegate<void(const GeometryContext& gctx,
                        const NavigationState& nState,
                        NavigationState::SurfaceCandidates& candidates),
                   ISurfaceCandidatesDelegate>;

/// Declare a Detctor Volume finding or switching delegate
///
/// @param gctx is the current geometry context
/// @param nState [in, out] is the navigation state to be updated
///
/// @return the new DetectorVolume into which one changes at this switch
using DetectorVolumeFinder =
    OwningDelegate<const DetectorVolume*(const GeometryContext& gctx,
                                         const NavigationState& nState),
                   IDetectorVolumeFinder>;

template <typename Derived, typename... Args>
inline static SurfaceCandidatesDelegate makeSurfaceCandidatesDelegate(
    Args... args) {
  SurfaceCandidatesDelegate delegate;
  delegate.template connect<&Derived::update>(
      std::make_unique<Derived>(std::forward<Args>(args)...));
  return delegate;
}

template <typename Derived, typename... Args>
inline static DetectorVolumeFinder makeDetectorVolumeFinder(Args... args) {
  DetectorVolumeFinder delegate;
  delegate.template connect<&Derived::find>(
      std::make_unique<Derived>(std::forward<Args>(args)...));
  return delegate;
}

}  // namespace Experimental
}  // namespace Acts
