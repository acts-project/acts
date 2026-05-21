// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// For whatever reason, this compilation unit does not compile
// with those assertions and GCC 13. For now just disable the
// flags in this case.
#if defined(_GLIBCXX_ASSERTIONS) && __GNUC__ == 13
#undef _GLIBCXX_ASSERTIONS
#endif

#include "ActsPlugins/Json/DetrayJsonHelper.hpp"

namespace Acts::DetrayJsonHelper {

std::size_t accelerationLink(std::span<const AxisDirection> casts) {
  // Default is `brute_force`
  using enum AxisDirection;
  std::size_t accLink = 0u;
  if (casts.size() == 2u) {
    if (casts[0u] == AxisX && casts[1u] == AxisY) {
      accLink = 1u;
    } else if (casts[0u] == AxisR && casts[1u] == AxisPhi) {
      accLink = 3u;
    } else if (casts[0u] == AxisZ && casts[1u] == AxisPhi) {
      accLink = 4u;
    } else if (casts[0u] == AxisZ && casts[1u] == AxisR) {
      accLink = 5u;
    }
  } else if (casts.size() == 3u) {
    if (casts[0u] == AxisX && casts[1u] == AxisY && casts[2u] == AxisZ) {
      accLink = 2u;
    } else if (casts[0u] == AxisZ && casts[1u] == AxisPhi &&
               casts[2u] == AxisR) {
      accLink = 5u;
    }
  }
  return accLink;
}
}  // namespace Acts::DetrayJsonHelper
