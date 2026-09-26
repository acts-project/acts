// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <utility>

namespace Acts {

class Surface;
class TrackingGeometry;
class TrackingVolume;

/// Type alias for surface material maps indexed by geometry identifier
using SurfaceMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
/// Type alias for volume material maps indexed by geometry identifier
using VolumeMaterialMaps =
    std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
/// A keyed assignment. The geometry ID records where it was mapped and is
/// diagnostic only when loading.
struct KeyedSurfaceMaterial {
  /// Geometry ID at mapping time, retained only for diagnostics.
  GeometryIdentifier geometryId;
  /// Material payload assigned by the stable key.
  std::shared_ptr<const ISurfaceMaterial> material;
};

/// Surface assignments with stable keys. Different keys may have the same
/// geometry ID recorded during mapping when maps from separate geometries are
/// combined.
using KeyedSurfaceMaterialMaps =
    std::map<std::string, KeyedSurfaceMaterial, std::less<>>;

/// Surface and volume material assignments.
/// Keyed assignments are stored separately and must never be matched by ID.
struct TrackingGeometryMaterial {
  /// Construct empty material assignments.
  TrackingGeometryMaterial() = default;

  /// Construct from material assignment maps.
  /// @param surfaces Unkeyed surface assignments
  /// @param volumes Volume assignments
  /// @param keyed Stable-key surface assignments
  TrackingGeometryMaterial(SurfaceMaterialMaps surfaces,
                           VolumeMaterialMaps volumes,
                           KeyedSurfaceMaterialMaps keyed = {})
      : surfaceMaterials(std::move(surfaces)),
        volumeMaterials(std::move(volumes)),
        keyedSurfaces(std::move(keyed)) {}

  /// Unkeyed surface assignments indexed by geometry ID.
  SurfaceMaterialMaps surfaceMaterials{};
  /// Volume assignments indexed by geometry ID.
  VolumeMaterialMaps volumeMaterials{};
  /// Surface assignments indexed by stable string key.
  KeyedSurfaceMaterialMaps keyedSurfaces{};

  /// Read the optional document description, ignored when applying material.
  /// @return Description, or std::nullopt if absent
  const std::optional<std::string>& description() const {
    return m_description;
  }

  /// Set or remove the document description. Empty descriptions are allowed.
  /// @param description Text, or std::nullopt to remove the description
  void setDescription(std::optional<std::string> description) {
    m_description = std::move(description);
  }

  /// Apply material to a completed geometry, checking surface identities and
  /// resolving all surface assignments before modifying the geometry.
  /// Keyed surfaces require a matching key; unused map entries are allowed.
  /// @throws std::invalid_argument if a keyed map is applied to Gen1 geometry
  /// @param geometry Geometry whose material is updated
  void apply(TrackingGeometry& geometry) const;

  /// Apply to selected surfaces with the same validation as a complete
  /// geometry.
  /// @param surfaces Surfaces whose material is updated; repeated pointers are allowed
  void apply(std::span<Surface* const> surfaces) const;

  /// Apply to a single surface. This cannot check uniqueness across surfaces.
  /// @param surface Surface whose material is updated
  void apply(Surface& surface) const;

  /// Apply volume material by geometry identifier.
  /// @param volume Volume whose material is updated
  void apply(TrackingVolume& volume) const;

 private:
  std::optional<std::string> m_description;
};

}  // namespace Acts
