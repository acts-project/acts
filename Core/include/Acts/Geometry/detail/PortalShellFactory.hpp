// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

namespace Acts {

class GeometryContext;
class PortalShellBase;
class TrackingVolume;

namespace detail {

/// Create the single-volume portal shell matching the bounds type of @p volume.
/// @param gctx The geometry context
/// @param volume The volume to create the shell for
/// @return The portal shell for the volume
/// @throws std::logic_error if the volume's bounds type has no shell
///         implementation
std::unique_ptr<PortalShellBase> makeSinglePortalShell(
    const GeometryContext& gctx, TrackingVolume& volume);

}  // namespace detail
}  // namespace Acts
