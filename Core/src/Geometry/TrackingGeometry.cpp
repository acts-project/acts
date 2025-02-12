// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>

namespace Acts {

class Gen1GeometryClosureVisitor : public TrackingGeometryMutableVisitor {
 public:
  Gen1GeometryClosureVisitor(const Logger& logger,
                             const GeometryIdentifierHook& hook)
      : m_logger(&logger), m_hook(&hook) {}

  const Logger& logger() const { return *m_logger; }

  void visitVolume(TrackingVolume& volume) override {
    // std::cout << "Volume: " << volume.volumeName() << std::endl;

    // Increment the volume ID for this volume
    m_volumeID = GeometryIdentifier().setVolume(m_volumeID.volume() + 1);
    // Reset boundary id for this volume
    m_iboundary = 0;
    // Reset layer id for this volume
    m_ilayer = 0;

    // std::cout << "volumeID: " << m_volumeID << std::endl;

    // assign the Volume ID to the volume itself
    volume.assignGeometryId(m_volumeID);
    ACTS_DEBUG("volumeID: " << m_volumeID << ", name: " << volume.volumeName());
    // insert the volume into the map
    m_volumesById[m_volumeID] = &volume;

    // assign the material if you have a decorator
    if (m_materialDecorator != nullptr) {
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
    // get the intersection solution
    auto& bSurface = boundary.surfaceRepresentation();
    // create the boundary surface id
    m_iboundary += 1;
    auto boundaryID = GeometryIdentifier(m_volumeID).setBoundary(m_iboundary);
    // std::cout << "boundaryID: " << boundaryID << std::endl;
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<RegularSurface*>(&bSurface));
    mutableBSurface.assignGeometryId(boundaryID);
    // Assign material if you have a decorator
    if (m_materialDecorator != nullptr) {
      m_materialDecorator->decorate(mutableBSurface);
    }
  }

  void visitLayer(Layer& layer) override {
    // create the layer identification
    m_ilayer += 1;
    auto layerID = GeometryIdentifier(m_volumeID).setLayer(m_ilayer);
    // now close the geometry
    layer.closeGeometry(m_materialDecorator, layerID, *m_hook, *m_logger);
  }

  void visitSurface(Surface& surface) override {
    if (surface.geometryId() == GeometryIdentifier{}) {
      throw std::invalid_argument("Surface has no geometry ID");
    }
    if (surface.geometryId().sensitive() != 0) {
      m_surfacesById[surface.geometryId()] = &surface;
    }
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

TrackingGeometry::TrackingGeometry(
    const MutableTrackingVolumePtr& highestVolume,
    const IMaterialDecorator* materialDecorator,
    const GeometryIdentifierHook& hook, const Logger& logger)
    : m_world(highestVolume) {
  Gen1GeometryClosureVisitor visitor{logger, hook};
  visitor.m_materialDecorator = materialDecorator;
  apply(visitor);

  m_volumesById = std::move(visitor.m_volumesById);
  m_surfacesById = std::move(visitor.m_surfacesById);
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

}  // namespace Acts
