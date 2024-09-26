// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NavigationPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

void TryAllPortalNavigationPolicy::updateState(
    const NavigationArguments& args) const {
  assert(m_volume != nullptr);

  for (const auto& portal : m_volume->portals()) {
    args.main.addPortalCandidate(portal);
  };
}
}  // namespace Acts
