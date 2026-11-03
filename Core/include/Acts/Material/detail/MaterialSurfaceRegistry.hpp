// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/TrackingGeometryMaterial.hpp"

#include <optional>
#include <span>

namespace Acts {
class Surface;
}

namespace Acts::detail {

/// Discover a stable assignment key on a designator placeholder.
const std::optional<std::string>& materialKey(const ISurfaceMaterial* material);

/// Validated identities of material surfaces in one completed geometry.
/// Repeated visits to the same surface are deduplicated. Distinct surfaces
/// must have distinct IDs and, where configured, distinct keys.
class MaterialSurfaceRegistry {
 public:
  /// Construct a registry for the participating material surfaces.
  /// @param surfaces Surfaces to validate, which must outlive the registry
  explicit MaterialSurfaceRegistry(std::span<const Surface* const> surfaces);

  /// Surfaces by geometry identifier
  const std::map<GeometryIdentifier, const Surface*>& surfaces() const {
    return m_surfaces;
  }

  /// Convert finalized ID-based mapping output to persistent assignments.
  /// @param materials Finalized mapping output
  TrackingGeometryMaterial materialMaps(SurfaceMaterialMaps materials) const;

 private:
  std::map<GeometryIdentifier, const Surface*> m_surfaces;
  std::map<GeometryIdentifier, std::optional<std::string>> m_keys;
};

}  // namespace Acts::detail
