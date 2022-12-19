// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Geometry/detail/NavigationStateUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

namespace Experimental {

struct DetectorShellBuilder {
  /// @brief
  /// @param shell
  /// @param gctx
  void expand(DetectorShell& shell, const GeometryContext& gctx) {}
};

}  // namespace Experimental

}  // namespace Acts
