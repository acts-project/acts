// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_map>

#include <detray/core/concepts.hpp>

namespace ActsPlugins {

/// @ingroup detray_plugin
/// @brief Two-way map between ACTS and detray surface identifiers
///
/// The detray surface identifiers are only assigned when the detray detector
/// is built, so the map is created from the built detector. The
/// @c DetrayPayloadConverter stores the ACTS geometry identifier of every
/// surface as its detray source link, which is what the map is built from.
///
/// Only sensitive surfaces are mapped: they have unique ACTS geometry
/// identifiers, while e.g. a single ACTS portal can be split into several
/// detray portals.
class DetraySurfaceIdMap {
 public:
  /// Type of a detray surface identifier (the detray barcode value)
  using DetrayId = std::uint64_t;

  /// Build the map from the sensitive surfaces of a detray detector
  /// @tparam detector_t the detray detector type
  /// @param detector the detray detector converted from ACTS
  /// @return the identifier map
  /// @throws std::invalid_argument if two sensitive surfaces share an ACTS
  ///         or a detray identifier
  template <detray::concepts::detector detector_t>
  static DetraySurfaceIdMap fromDetector(const detector_t& detector) {
    DetraySurfaceIdMap map;
    for (const auto& surface : detector.surfaces()) {
      if (!surface.is_sensitive()) {
        continue;
      }
      map.insert(surface.identifier().value(),
                 Acts::GeometryIdentifier(surface.source));
    }
    return map;
  }

  /// Add one surface to the map
  /// @param detrayId the detray identifier of the surface
  /// @param id the ACTS geometry identifier of the surface
  /// @throws std::invalid_argument if either identifier is already mapped
  void insert(DetrayId detrayId, Acts::GeometryIdentifier id);

  /// detray identifier of an ACTS surface
  /// @param id the ACTS geometry identifier
  /// @return the detray identifier, or nothing if the surface is not mapped
  std::optional<DetrayId> toDetray(Acts::GeometryIdentifier id) const;

  /// ACTS geometry identifier of a detray surface
  /// @param detrayId the detray identifier
  /// @return the ACTS geometry identifier, or nothing if the surface is not
  ///         mapped
  std::optional<Acts::GeometryIdentifier> fromDetray(DetrayId detrayId) const;

  /// @return the number of mapped surfaces
  std::size_t size() const { return m_fromDetray.size(); }

  /// @return true if no surface is mapped
  bool empty() const { return m_fromDetray.empty(); }

 private:
  std::unordered_map<DetrayId, Acts::GeometryIdentifier> m_fromDetray;
  std::unordered_map<Acts::GeometryIdentifier, DetrayId> m_toDetray;
};

}  // namespace ActsPlugins
