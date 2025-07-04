// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

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
    if (volume.volumeMaterial() == nullptr &&
        volume.motherVolume() != nullptr &&
        volume.motherVolume()->volumeMaterial() != nullptr) {
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
                << surface.toStream(GeometryContext()) << std::endl;
      throw std::invalid_argument("Surface has no geometry ID");
    }

    checkIdentifier(surface, "surface");

    //@TODO: Why not use all of them?
    if (surface.geometryId().sensitive() != 0) {
      m_surfacesById.emplace(surface.geometryId(), &surface);
    }
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
    checkIdentifier(boundary.surfaceRepresentation(), "boundary surface");
  }

  void visitPortal(const Portal& portal) override {
    checkIdentifier(portal.surface(), "portal");
  }

  std::unordered_map<GeometryIdentifier, const TrackingVolume*> m_volumesById{};
  std::unordered_map<GeometryIdentifier, const Surface*> m_surfacesById{};

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

  ACTS_DEBUG("TrackingGeometry created with "
             << m_volumesById.size() << " volumes and " << m_surfacesById.size()
             << " sensitive surfaces");

  m_volumesById.rehash(0);
  m_surfacesById.rehash(0);
}

TrackingGeometry::~TrackingGeometry() = default;

const TrackingVolume* TrackingGeometry::lowestTrackingVolume(
    const GeometryContext& gctx, const Vector3& gp) const {
  return m_world->lowestTrackingVolume(gctx, gp, s_onSurfaceTolerance);
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

const Surface* TrackingGeometry::findSurface(GeometryIdentifier id) const {
  auto srf = m_surfacesById.find(id);
  if (srf == m_surfacesById.end()) {
    return nullptr;
  }
  return srf->second;
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
