// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/GeometryIdGenerators.hpp"

#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <sstream>

void Acts::Experimental::detail::VolumeCounter::generateIds(
    DetectorVolume& volume) {
  // The geometry Id has already been set, hence ignore
  if (volume.geometryId().volume() != 0u) {
    return;
  }

  GeometryIdentifier volumeID;
  volumeID.setVolume(++volumeCounter);
  volume.assignGeometryId(volumeID);

  if (processSurfacesPortals) {
    // Set to all contained surfaces
    size_t sensitiveCounter = 0;
    size_t passiveCounter = 0;
    for (auto sPtr : volume.surfacePtrs()) {
      auto surfaceID = sPtr->geometryId();
      if (sPtr->associatedDetectorElement() != nullptr and
          surfaceID.sensitive() == 0u) {
        sPtr->assignGeometryId(
            (GeometryIdentifier(volumeID).setSensitive(++sensitiveCounter)));
      } else if (surfaceID.passive() == 0u and
                 sPtr->associatedDetectorElement() == nullptr)
        sPtr->assignGeometryId(
            (GeometryIdentifier(volumeID).setPassive(++passiveCounter)));
    }

    // Set to the portal if not yet claimed by another volume
    size_t portalCounter = 0;
    for (auto pPtr : volume.portalPtrs()) {
      if (pPtr->surface().geometryId().boundary() == 0u) {
        // Claim the portal surface
        pPtr->surface().assignGeometryId(
            (GeometryIdentifier(volumeID).setBoundary(++portalCounter)));
      }
    }
  }

  // Recursively step down down into sub volumes if present
  for (auto vPtr : volume.volumePtrs()) {
    generateIds(*vPtr);
  }
}

void Acts::Experimental::detail::LayerCounter::generateIds(
    DetectorVolume& volume) {
  // The geometry Id has already been set, hence ignore
  if (volume.geometryId().volume() != 0u) {
    return;
  }
  volume.assignGeometryId(GeometryIdentifier(baseID).setLayer(++layerCounter));
}

void Acts::Experimental::detail::PortalCounter::generateIds(
    DetectorVolume& volume) {
  size_t portalCounter = 0;
  for (auto pPtr : volume.portalPtrs()) {
    if (pPtr->surface().geometryId().boundary() == 0u) {
      // Claim the portal surface
      pPtr->surface().assignGeometryId(
          GeometryIdentifier(volume.geometryId()).setBoundary(++portalCounter));
    }
  }
}

void Acts::Experimental::detail::SensitiveCounter::generateIds(
    DetectorVolume& volume) {
  size_t sensitiveCounter = 0;
  for (auto& sPtr : volume.surfacePtrs()) {
    auto surfaceID = sPtr->geometryId();
    // Set only if the surface is indeed sensitive and if the sensitive flag
    // is not yet set
    if (sPtr->associatedDetectorElement() != nullptr and
        surfaceID.sensitive() == 0u) {
      sPtr->assignGeometryId(GeometryIdentifier(volume.geometryId())
                                 .setSensitive(++sensitiveCounter));
    }
  }
}

void Acts::Experimental::detail::PassiveCounter::generateIds(
    DetectorVolume& volume) {
  size_t passiveCounter = 0;
  for (auto sPtr : volume.surfacePtrs()) {
    auto surfaceID = sPtr->geometryId();
    // Set only if the surface is indeed passive
    if (surfaceID.passive() == 0u and
        sPtr->associatedDetectorElement() == nullptr) {
      sPtr->assignGeometryId(
          GeometryIdentifier(volume.geometryId()).setPassive(++passiveCounter));
    }
  }
}

void Acts::Experimental::detail::DuplicateIdChecker::generateIds(
    DetectorVolume& volume) {
  // Check the geometry id values
  auto check = [&](const GeometryIdentifier& geoID) -> void {
    if (geoID.value() == 0u) {
      // This needs be found by the unset Id checker
      return;
    }

    if (m_geometryIdMap.find(geoID) != m_geometryIdMap.end()) {
      std::stringstream ss;
      ss << "DuplicateIdChecker: duplicate GeometryIdentifier ";
      ss << geoID << " in volume " << volume.geometryId();
      ss << " found.";
      throw std::runtime_error(ss.str());
    }
    m_geometryIdMap[geoID] = true;
  };
  // Check the volume itself
  check(volume.geometryId());
  m_processedVolumes[&volume] = true;
  // Check contained surfaces
  for (auto s : volume.surfaces()) {
    check(s->geometryId());
  }
  // Check contained portals, but keep track of shared ones
  for (auto p : volume.portals()) {
    if (m_processedPortals.find(p) == m_processedPortals.end()) {
      check(p->surface().geometryId());
      m_processedPortals[p] = true;
    }
  }
  // Check contained volumes, allow for global appearance in detector listx
  for (auto v : volume.volumePtrs()) {
    if (m_processedVolumes.find(v.get()) == m_processedVolumes.end()) {
      generateIds(*v);
      m_processedVolumes[v.get()] = true;
    }
  }
}

void Acts::Experimental::detail::UnsetIdChecker::generateIds(
    DetectorVolume& volume) {
  auto check = [&](const GeometryIdentifier& geoID) -> void {
    if (geoID.value() == 0u) {
      std::stringstream ss;
      ss << "UnsetIdChecker: unset GeometryIdentifier ";
      ss << geoID << " in volume " << volume.geometryId();
      ss << " found.";
      throw std::runtime_error(ss.str());
    }
  };
  // Check the volume itself
  check(volume.geometryId());
  // Check the surfaces
  for (auto s : volume.surfaces()) {
    check(s->geometryId());
  }
  // Check the portals
  for (auto p : volume.portals()) {
    check(p->surface().geometryId());
  }
  // Recursively step down
  for (auto v : volume.volumePtrs()) {
    generateIds(*v);
  }
}