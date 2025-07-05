// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayFwd.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"

#include <memory>

#include <detray/io/frontend/payloads.hpp>

namespace Acts {

std::unique_ptr<DetraySurfaceGrid> MultiNavigationPolicy::toDetrayPayload()
    const {
  // Only ONE of the child policies should return a non-nullptr payload
  for (const auto& policy : m_policyPtrs) {
    auto payload = policy->toDetrayPayload();
    if (payload) {
      return payload;
    }
  }
  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid>
Experimental::MultiLayerNavigationPolicy::toDetrayPayload() const {
  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid>
SurfaceArrayNavigationPolicy::toDetrayPayload() const {
  // @TODO: Implement this conversion!
  return nullptr;
}

std::unique_ptr<DetraySurfaceGrid> TryAllNavigationPolicy::toDetrayPayload()
    const {
  return nullptr;
}

}  // namespace Acts
