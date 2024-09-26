// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

namespace Experimental::Gen3Geometry {
// @NOTE: This type is PRELIMINARY! It does not represent the final
//        implementation of the updated navigation!
struct NavigationState {
  const TrackingVolume* currentVolume = nullptr;
  NavigationStream main;
};

}  // namespace Experimental::Gen3Geometry

class INavigationPolicy {
 public:
  // @NOTE: This interface is PRELIMINARY! It is subject to change!

  using DelegateType =
      Delegate<void(Experimental::Gen3Geometry::NavigationState&)>;

  static void noopUpdate(
      Experimental::Gen3Geometry::NavigationState& /*state*/) {}

  virtual ~INavigationPolicy() = default;

  virtual void connect(DelegateType& delegate) const = 0;
};

class TryAllPortalNavigationPolicy final : public INavigationPolicy {
 public:
  void connect(DelegateType& delegate) const override {
    delegate.connect<&TryAllPortalNavigationPolicy::updateState>(this);
  }

  TryAllPortalNavigationPolicy(const TrackingVolume& volume)
      : m_volume(&volume) {}

 private:
  // @NOTE: This implementation is PRELIMINARY! It is subject to change!
  void updateState(Experimental::Gen3Geometry::NavigationState& state) const;

  const TrackingVolume* m_volume;
};

}  // namespace Acts
