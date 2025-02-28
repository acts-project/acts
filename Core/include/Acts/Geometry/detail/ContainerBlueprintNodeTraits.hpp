// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeStack.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"

namespace Acts::Experimental {

template <class Bounds>
struct ContainerBlueprintNodeTraits;

template <>
struct ContainerBlueprintNodeTraits<CylinderVolumeBounds> {
  using VolumeStack = CylinderVolumeStack;

  using BaseShell = CylinderPortalShell;
  using SingleShell = SingleCylinderPortalShell;
  using ShellStack = CylinderStackPortalShell;
};

template <>
struct ContainerBlueprintNodeTraits<CuboidVolumeBounds> {
  using VolumeStack = CuboidVolumeStack;

  using BaseShell = CuboidPortalShell;
  using SingleShell = SingleCuboidPortalShell;
  using ShellStack = CuboidStackPortalShell;
};

}  // namespace Acts::Experimental
