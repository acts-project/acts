// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>
#include <format>
#include <optional>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <detray/core/concepts.hpp>
#include <detray/geometry/identifier.hpp>

namespace ActsPlugins {

/// @ingroup detray_plugin
/// @brief Two-way map between ACTS and detray surface identifiers
///
/// The detray surface identifiers are only assigned when the detray detector
/// is built, so the map is created from the built detector. The
/// @c DetrayPayloadConverter stores the ACTS geometry identifier of every
/// surface as its detray source link, which is what the map is built from.
///
/// Besides the identifiers, the map also gives direct access to the ACTS
/// surface of every mapped detray surface. Device code cannot hold pointers,
/// so the detray identifier stays the key that is passed around, e.g. as the
/// surface link of traccc track parameters.
///
/// Only sensitive surfaces are mapped: they have unique ACTS geometry
/// identifiers, while e.g. a single ACTS portal can be split into several
/// detray portals.
///
/// The map cannot be changed after construction. Only the surface transforms
/// change with alignment, and those are handled through the geometry context.
///
/// @note The map holds non-owning pointers to the ACTS surfaces, so the
///       tracking geometry has to outlive it.
class DetraySurfaceIdMap {
 public:
  /// One mapped surface: its detray identifier and the ACTS surface
  using Entry = std::pair<detray::geometry::identifier, const Acts::Surface*>;

  /// Create an empty map
  DetraySurfaceIdMap() = default;

  /// Build the map from a list of surfaces
  /// @param entries the detray identifier and ACTS surface of every surface
  /// @throws std::invalid_argument if a surface is null, or if two entries
  ///         share an ACTS or a detray identifier
  explicit DetraySurfaceIdMap(std::span<const Entry> entries);

  /// Build the map from the sensitive surfaces of a detray detector
  /// @tparam detector_t the detray detector type
  /// @param detector the detray detector converted from ACTS
  /// @param trackingGeometry the ACTS tracking geometry @p detector was
  ///        converted from
  /// @throws std::invalid_argument if a sensitive detray surface has no ACTS
  ///         surface in @p trackingGeometry, or if two sensitive surfaces
  ///         share an ACTS or a detray identifier
  template <detray::concepts::detector detector_t>
  DetraySurfaceIdMap(const detector_t& detector,
                     const Acts::TrackingGeometry& trackingGeometry)
      : DetraySurfaceIdMap(collectEntries(detector, trackingGeometry)) {}

  /// detray identifier of an ACTS surface
  /// @param id the ACTS geometry identifier
  /// @return the detray identifier, or nothing if the surface is not mapped
  std::optional<detray::geometry::identifier> toDetray(
      Acts::GeometryIdentifier id) const;

  /// ACTS geometry identifier of a detray surface
  /// @param detrayId the detray identifier
  /// @return the ACTS geometry identifier, or nothing if the surface is not
  ///         mapped
  std::optional<Acts::GeometryIdentifier> fromDetray(
      detray::geometry::identifier detrayId) const;

  /// ACTS surface of a detray surface
  /// @param detrayId the detray identifier
  /// @return the ACTS surface, or nullptr if the surface is not mapped
  const Acts::Surface* surface(detray::geometry::identifier detrayId) const;

  /// @return the number of mapped surfaces
  std::size_t size() const { return m_fromDetray.size(); }

  /// @return true if no surface is mapped
  bool empty() const { return m_fromDetray.empty(); }

 private:
  template <detray::concepts::detector detector_t>
  static std::vector<Entry> collectEntries(
      const detector_t& detector,
      const Acts::TrackingGeometry& trackingGeometry) {
    std::vector<Entry> entries;
    for (const auto& surface : detector.surfaces()) {
      if (!surface.is_sensitive()) {
        continue;
      }
      const Acts::GeometryIdentifier id(surface.source);
      const Acts::Surface* actsSurface = trackingGeometry.findSurface(id);
      if (actsSurface == nullptr) {
        throw std::invalid_argument(
            std::format("DetraySurfaceIdMap: no ACTS surface with id {} for "
                        "detray surface {}",
                        id, surface.identifier().value()));
      }
      entries.emplace_back(surface.identifier(), actsSurface);
    }
    return entries;
  }

  std::unordered_map<detray::geometry::identifier, const Acts::Surface*>
      m_fromDetray;
  std::unordered_map<Acts::GeometryIdentifier, detray::geometry::identifier>
      m_toDetray;
};

}  // namespace ActsPlugins
