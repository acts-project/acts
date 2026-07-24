// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cassert>
#include <cstddef>

namespace Acts {

class Gen1GeometryClosureVisitor : public TrackingGeometryMutableVisitor {
 public:
  Gen1GeometryClosureVisitor(const Logger& logger,
                             const IMaterialDecorator* materialDecorator,
                             const GeometryIdentifierHook& hook)
      : m_logger(&logger),
        m_materialDecorator(materialDecorator),
        m_hook(&hook) {
    ACTS_VERBOSE("Creating Gen1GeometryClosureVisitor");
  }

  const Logger& logger() const { return *m_logger; }

  void visitVolume(TrackingVolume& volume) override {
    ACTS_DEBUG("Volume: " << volume.volumeName());

    // Increment the volume ID for this volume
    m_volumeID = GeometryIdentifier().withVolume(m_volumeID.volume() + 1);
    // Reset boundary id for this volume
    m_iboundary = 0;
    // Reset layer id for this volume
    m_ilayer = 0;

    // assign the Volume ID to the volume itself
    ACTS_VERBOSE("~> volumeID: " << m_volumeID);
    volume.assignGeometryId(m_volumeID);

    // assign the material if you have a decorator
    if (m_materialDecorator != nullptr) {
      ACTS_VERBOSE("Decorating volume " << volume.volumeName()
                                        << " with material");
      m_materialDecorator->decorate(volume);
    }
    if (!volume.hasMaterial() && volume.motherVolume() != nullptr &&
        volume.motherVolume()->hasMaterial()) {
      auto protoMaterial = dynamic_cast<const ProtoVolumeMaterial*>(
          volume.motherVolume()->volumeMaterial());
      if (protoMaterial == nullptr) {
        volume.assignVolumeMaterial(volume.motherVolume()->volumeMaterialPtr());
      }
    }
  }

  void visitBoundarySurface(
      BoundarySurfaceT<TrackingVolume>& boundary) override {
    ACTS_DEBUG("BoundarySurface: " << boundary.surfaceRepresentation().name());
    // get the intersection solution
    auto& bSurface = boundary.surfaceRepresentation();
    // create the boundary surface id
    m_iboundary += 1;
    auto boundaryID = GeometryIdentifier(m_volumeID).withBoundary(m_iboundary);
    ACTS_VERBOSE("~> boundaryID: " << boundaryID);
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<RegularSurface*>(&bSurface));

    // assign the boundary ID to the surface
    ACTS_VERBOSE("~> assigning boundaryID: " << boundaryID);
    mutableBSurface.assignGeometryId(boundaryID);

    // Assign material if you have a decorator
    if (m_materialDecorator != nullptr) {
      ACTS_VERBOSE("Decorating boundary surface " << bSurface.name()
                                                  << " with material");
      m_materialDecorator->decorate(mutableBSurface);
    }
  }

  void visitLayer(Layer& layer) override {
    ACTS_DEBUG("Close Layer");
    // create the layer identification
    m_ilayer += 1;
    auto layerID = GeometryIdentifier(m_volumeID).withLayer(m_ilayer);
    ACTS_VERBOSE("~> layerID: " << layerID);

    // now close the geometry
    layer.closeGeometry(m_materialDecorator, layerID, *m_hook, *m_logger);
  }

  const Logger* m_logger;
  GeometryIdentifier m_volumeID;
  GeometryIdentifier::Value m_iboundary = 0;
  GeometryIdentifier::Value m_ilayer = 0;
  const IMaterialDecorator* m_materialDecorator = nullptr;
  const GeometryIdentifierHook* m_hook = nullptr;

  std::unordered_map<GeometryIdentifier, const TrackingVolume*> m_volumesById{};
  std::unordered_map<GeometryIdentifier, const Surface*> m_surfacesById{};
};

namespace {
class GeometryIdMapVisitor : public TrackingGeometryVisitor {
 private:
  void checkIdentifier(const GeometryObject& obj, std::string_view type) {
    if (obj.geometryId() == GeometryIdentifier{}) {
      std::stringstream ss;
      ss << "Encountered " << type << " with no geometry ID";
      throw std::invalid_argument(ss.str());
    }

    ACTS_VERBOSE("Checking identifier for " << type << ": "
                                            << obj.geometryId());

    auto [it, inserted] = m_objectsById.emplace(obj.geometryId(), &obj);

    if (!inserted && it->second != &obj) {
      std::stringstream ss;
      ss << "Duplicate " << type << " ID: " << obj.geometryId() << ": & "
         << it->second << " != " << &obj;
      if (const auto* other = dynamic_cast<const TrackingVolume*>(it->second);
          other != nullptr) {
        ss << " (" << other->volumeName() << ")";
      }
      ACTS_ERROR(ss.str());
      throw std::invalid_argument(ss.str());
    } else {
      ACTS_VERBOSE("Inserted " << type << " ID: " << obj.geometryId()
                               << " pointing at " << &obj);
    }
  }

  const Logger& logger() const { return m_logger; }
  const Logger& m_logger;

 public:
  explicit GeometryIdMapVisitor(const Logger& logger) : m_logger(logger) {}

  void visitVolume(const TrackingVolume& volume) override {
    std::string label = "volume(" + volume.volumeName() + ")";
    checkIdentifier(volume, label);

    m_volumesById.emplace(volume.geometryId(), &volume);
  }

  void visitSurface(const Surface& surface) override {
    if (surface.geometryId() == GeometryIdentifier{}) {
      std::cout << "Surface has no geometry ID: "
                << surface.toStream(
                       GeometryContext::dangerouslyDefaultConstruct())
                << std::endl;
      throw std::invalid_argument("Surface has no geometry ID");
    }

    checkIdentifier(surface, "surface");

    m_surfacesById.emplace(surface.geometryId(), &surface);
  }

  void visitLayer(const Layer& layer) override {
    // Layers ARE also GeometryObjects and have IDs.
    // Let's check that the layer has the same ID as it's surface
    // representation. Uniqueness of the surface IDs is checked in the surface
    if (layer.geometryId() != layer.surfaceRepresentation().geometryId()) {
      ACTS_ERROR("Layer ID mismatch: "
                 << layer.geometryId()
                 << " != " << layer.surfaceRepresentation().geometryId());
      throw std::invalid_argument("Layer ID mismatch");
    }
  }

  void visitBoundarySurface(
      const BoundarySurfaceT<TrackingVolume>& boundary) override {
    const auto& surface = boundary.surfaceRepresentation();
    checkIdentifier(surface, "boundary surface");
    m_surfacesById.emplace(surface.geometryId(), &surface);
  }

  void visitPortal(const Portal& portal) override {
    const auto& surface = portal.surface();
    checkIdentifier(surface, "portal");
    m_surfacesById.emplace(surface.geometryId(), &surface);

    for (const auto& tag : portal.tags()) {
      auto [it, inserted] = m_portalsByTag.try_emplace(tag, &portal);
      // A fused/merged portal is shared between volumes, so it is visited once
      // per owning volume slot. Re-inserting the same tag for the *same* portal
      // is fine; a different portal claiming the same tag is a collision.
      if (!inserted && it->second != &portal) {
        std::stringstream ss;
        ss << "Duplicate portal tag: " << tag;
        ACTS_ERROR(ss.str());
        throw std::invalid_argument(ss.str());
      }
    }
  }

  std::unordered_map<GeometryIdentifier, const TrackingVolume*> m_volumesById{};
  std::unordered_map<GeometryIdentifier, const Surface*> m_surfacesById{};
  detail::PortalTagMap m_portalsByTag{};

  std::unordered_map<GeometryIdentifier, const GeometryObject*> m_objectsById{};
};

}  // namespace
TrackingGeometry::TrackingGeometry(
    const MutableTrackingVolumePtr& highestVolume,
    const IMaterialDecorator* materialDecorator,
    const GeometryIdentifierHook& hook, const Logger& logger, bool close)
    : m_world(highestVolume) {
  if (close) {
    ACTS_DEBUG("Closing tracking geometry with Gen1 assignment");
    Gen1GeometryClosureVisitor visitor{logger, materialDecorator, hook};
    apply(visitor);
  }

  GeometryIdMapVisitor mapVisitor{logger};
  apply(mapVisitor);
  m_volumesById = std::move(mapVisitor.m_volumesById);
  m_surfacesById = std::move(mapVisitor.m_surfacesById);
  m_portalsByTag = std::move(mapVisitor.m_portalsByTag);

  ACTS_DEBUG("TrackingGeometry created with "
             << m_volumesById.size() << " volumes and " << m_surfacesById.size()
             << " surfaces");

  m_volumesById.rehash(0);
  m_surfacesById.rehash(0);
  m_portalsByTag.rehash(0);
}

TrackingGeometry::~TrackingGeometry() = default;

namespace {

/// Resolve the volume a track on @p surface at @p position is entering,
/// assuming it is currently associated with @p volume. If the surface is one
/// of the boundary surfaces (Gen1) or portal surfaces (Gen3) of the volume,
/// the volume on the far side along @p direction is returned, which can be
/// `nullptr` if there is no volume in that direction (end of world).
/// Otherwise the surface does not act as a boundary here and @p volume
/// itself is returned.
///
/// The position is assumed to be on @p surface. A position inside the
/// volume and on the plane of a matched boundary is within that boundary's
/// bounds up to the lookup tolerance, since portals and boundary surfaces
/// cover the volume faces. The exact bounds check can still fail for
/// positions grazing a volume edge within the tolerance; in that case the
/// boundary is not crossed here and the volume itself is returned.
const TrackingVolume* resolveVolumeThroughBoundary(const GeometryContext& gctx,
                                                   const TrackingVolume& volume,
                                                   const Vector3& position,
                                                   double tolerance,
                                                   const Vector3& direction,
                                                   const Surface& surface) {
  auto isOnBoundary = [&]() {
    return surface.isOnSurface(gctx, position, direction,
                               BoundaryTolerance::None(), tolerance);
  };

  for (const Portal& portal : volume.portals()) {
    if (&portal.surface() == &surface) {
      if (!isOnBoundary()) {
        return &volume;
      }
      auto resolved = portal.resolveVolume(gctx, position, direction);
      if (!resolved.ok()) {
        return &volume;
      }
      return resolved.value();
    }
  }

  for (const auto& boundary : volume.boundarySurfaces()) {
    if (&boundary->surfaceRepresentation() == &surface) {
      if (!isOnBoundary()) {
        return &volume;
      }
      return boundary->attachedVolume(gctx, position, direction);
    }
  }

  return &volume;
}

}  // namespace

const TrackingVolume* TrackingGeometry::lowestTrackingVolume(
    const GeometryContext& gctx, const Vector3& gp, double tolerance,
    const std::optional<Vector3>& direction,
    const Surface* associatedSurface) const {
  const TrackingVolume* volume =
      m_world->lowestTrackingVolume(gctx, gp, tolerance);
  if (volume != nullptr && direction.has_value() &&
      associatedSurface != nullptr) {
    assert(associatedSurface->isOnSurface(gctx, gp, *direction,
                                          BoundaryTolerance::Infinite(),
                                          tolerance) &&
           "The associated surface must contain the position");
    volume = resolveVolumeThroughBoundary(gctx, *volume, gp, tolerance,
                                          *direction, *associatedSurface);
  }
  return volume;
}

const TrackingVolume* TrackingGeometry::highestTrackingVolume() const {
  return m_world.get();
}

TrackingVolume* TrackingGeometry::highestTrackingVolume() {
  return m_world.get();
}

std::shared_ptr<const TrackingVolume>
TrackingGeometry::highestTrackingVolumePtr() const {
  return m_world;
}

const Layer* TrackingGeometry::associatedLayer(const GeometryContext& gctx,
                                               const Vector3& gp) const {
  const TrackingVolume* lowestVol = lowestTrackingVolume(gctx, gp);
  if (lowestVol == nullptr) {
    return nullptr;
  }
  return lowestVol->associatedLayer(gctx, gp);
}

const TrackingVolume* TrackingGeometry::findVolume(
    GeometryIdentifier id) const {
  auto vol = m_volumesById.find(id);
  if (vol == m_volumesById.end()) {
    return nullptr;
  }
  return vol->second;
}

const TrackingVolume* TrackingGeometry::findVolumeByName(
    std::string_view name) const {
  const TrackingVolume* found = nullptr;
  apply([&](const TrackingVolume& volume) {
    if (found == nullptr && volume.volumeName() == name) {
      found = &volume;
    }
  });
  return found;
}

const Surface* TrackingGeometry::findSurface(GeometryIdentifier id) const {
  auto srf = m_surfacesById.find(id);
  if (srf == m_surfacesById.end()) {
    return nullptr;
  }
  return srf->second;
}

const Portal* TrackingGeometry::findPortal(std::string_view tag) const {
  auto it = m_portalsByTag.find(tag);
  if (it == m_portalsByTag.end()) {
    return nullptr;
  }
  return it->second;
}

const std::unordered_map<GeometryIdentifier, const Surface*>&
TrackingGeometry::geoIdSurfaceMap() const {
  return m_surfacesById;
}

void TrackingGeometry::visualize(IVisualization3D& helper,
                                 const GeometryContext& gctx,
                                 const ViewConfig& viewConfig,
                                 const ViewConfig& portalViewConfig,
                                 const ViewConfig& sensitiveViewConfig) const {
  highestTrackingVolume()->visualize(helper, gctx, viewConfig, portalViewConfig,
                                     sensitiveViewConfig);
}

void TrackingGeometry::apply(TrackingGeometryVisitor& visitor) const {
  highestTrackingVolume()->apply(visitor);
}

void TrackingGeometry::apply(TrackingGeometryMutableVisitor& visitor) {
  highestTrackingVolume()->apply(visitor);
}

TrackingGeometry::GeometryVersion TrackingGeometry::geometryVersion() const {
  if (highestTrackingVolume()->portals().empty()) {
    return GeometryVersion::Gen1;
  } else {
    return GeometryVersion::Gen3;
  }
}

}  // namespace Acts
