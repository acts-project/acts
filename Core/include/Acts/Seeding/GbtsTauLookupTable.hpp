// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts::Experimental {

/// Accepted |cot(theta)| range for one cluster width bin. A cluster near the
/// module edge may be shortened, hence the second pair.
struct GbtsTauBounds final {
  /// Minimum accepted |cot(theta)|.
  float minTau{};
  /// Maximum accepted |cot(theta)|. Negative means undertrained: do not cut.
  float maxTau{};
  /// Minimum accepted |cot(theta)| near the module edge.
  float minTauNearEdge{};
  /// Maximum accepted |cot(theta)| near the module edge. Negative as above.
  float maxTauNearEdge{};
};

/// Tau bounds per cluster width, one entry per 0.05 mm, indexed not searched.
using GbtsTauLookupTable = std::vector<GbtsTauBounds>;

}  // namespace Acts::Experimental
