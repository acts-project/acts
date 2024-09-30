// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Navigation/INavigationPolicy.hpp"

namespace Acts {

class TrackingVolume;

// @NOTE: This implementation is PRELIMINARY! It is subject to change!

class TryAllPortalNavigationPolicy final : public INavigationPolicy {
 public:
  explicit TryAllPortalNavigationPolicy(const TrackingVolume& volume);

  void updateState(const NavigationArguments& args) const;

  void connect(NavigationDelegate& delegate) const override;

 private:
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllPortalNavigationPolicy>);

class TryAllSurfaceNavigationPolicy final : public INavigationPolicy {
 public:
  explicit TryAllSurfaceNavigationPolicy(const TrackingVolume& volume);

  void updateState(const NavigationArguments& args) const;

  void connect(NavigationDelegate& delegate) const override;

 private:
  const TrackingVolume* m_volume;
};

static_assert(NavigationPolicyConcept<TryAllSurfaceNavigationPolicy>);

}  // namespace Acts
