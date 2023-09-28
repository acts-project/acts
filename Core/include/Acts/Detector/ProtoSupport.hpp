// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <optional>
#include <stdexcept>
#include <vector>

namespace Acts {

namespace Experimental {
/// @brief Support parameter definitions
struct ProtoSupport {
  /// Define whether you want to build support structures
  std::array<ActsScalar, 5u> values = {};
  /// The surface type to be built
  Surface::SurfaceType type = Surface::SurfaceType::Other;
  /// Define in which values the support should be constrained
  std::vector<BinningValue> constraints = s_binningValues;
  /// Potential splits into planar approximations
  unsigned int splits = 1u;
  /// The (optional) layer transform
  std::optional<Transform3> transform = std::nullopt;
  /// Alternatively - the support surface can already be provided
  std::shared_ptr<Surface> surface = nullptr;
  /// Indicate if the support surface should always be addressed in navigation
  bool assignToAll = false;
};

}  // namespace Experimental
}  // namespace Acts
