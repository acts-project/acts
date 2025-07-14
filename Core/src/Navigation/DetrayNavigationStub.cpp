// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayExceptions.hpp"
#include "Acts/Geometry/DetrayFwd.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"

namespace Acts {

#define STUB_METHOD(type)                                   \
  std::unique_ptr<DetraySurfaceGrid> type::toDetrayPayload( \
      const SurfaceLookupFunction& surfaceLookup) const {   \
    throw DetrayNotAvailableException();                    \
  }

// In STUB mode: all navigation related methods throw an exception to indicate
// that Detray is not available.

// clang-format off

STUB_METHOD(MultiNavigationPolicy)
STUB_METHOD(SurfaceArrayNavigationPolicy)
STUB_METHOD(Experimental::MultiLayerNavigationPolicy)
STUB_METHOD(TryAllNavigationPolicy)

// clang-format on

}  // namespace Acts
