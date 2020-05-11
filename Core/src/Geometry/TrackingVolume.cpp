// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <functional>
#include <utility>

#include "Acts/Geometry/GlueVolumesDescriptor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"

Acts::TrackingVolume::TrackingVolume()
    : Volume(),
      m_volumeMaterial(nullptr),
      m_boundarySurfaces(),
      m_confinedLayers(nullptr),
      m_confinedVolumes(nullptr),
      m_name("undefined") {}

Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volbounds,
    const std::shared_ptr<const TrackingVolumeArray>& containedVolumeArray,
    const std::string& volumeName)
    : Volume(std::move(htrans), std::move(volbounds)),
      m_volumeMaterial(nullptr),
      m_boundarySurfaces(),
      m_confinedLayers(nullptr),
      m_confinedVolumes(containedVolumeArray),
      m_name(volumeName) {
  createBoundarySurfaces();
  interlinkLayers();
}

// constructor for arguments
Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volumeBounds,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    std::unique_ptr<const LayerArray> staticLayerArray,
    std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
    MutableTrackingVolumeVector denseVolumeVector,
    const std::string& volumeName)
    : Volume(std::move(htrans), std::move(volumeBounds)),
      m_volumeMaterial(std::move(volumeMaterial)),
      m_confinedLayers(std::move(staticLayerArray)),
      m_confinedVolumes(std::move(containedVolumeArray)),
      m_confinedDenseVolumes({}),
      m_name(volumeName) {
  createBoundarySurfaces();
  interlinkLayers();
  connectDenseBoundarySurfaces(denseVolumeVector);
}

// constructor for arguments
Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volbounds,
    std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
    std::vector<std::unique_ptr<const Volume>> descendants,
    const Volume::BoundingBox* top,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    const std::string& volumeName)
    : Volume(std::move(htrans), std::move(volbounds)),
      m_volumeMaterial(std::move(volumeMaterial)),
      m_name(volumeName),
      m_descendantVolumes(std::move(descendants)),
      m_bvhTop(top) {
  createBoundarySurfaces();
  // we take a copy of the unique box pointers, but we want to
  // store them as consts.
  for (auto& uptr : boxStore) {
    m_boundingBoxes.push_back(
        std::unique_ptr<Volume::BoundingBox>(uptr.release()));
  }
}

Acts::TrackingVolume::~TrackingVolume() {
  delete m_glueVolumeDescriptor;
}

const Acts::TrackingVolume* Acts::TrackingVolume::lowestTrackingVolume(
    const GeometryContext& /*gctx*/, const Vector3D& position,
    const double tol) const {
  // confined static volumes - highest hierarchy
  if (m_confinedVolumes) {
    return (m_confinedVolumes->object(position).get());
  }

  // search for dense volumes
  if (!m_confinedDenseVolumes.empty())
    for (auto& denseVolume : m_confinedDenseVolumes)
      if (denseVolume->inside(position, tol))
        return denseVolume.get();

  // there is no lower sub structure
  return this;
}

const Acts::TrackingVolumeBoundaries& Acts::TrackingVolume::boundarySurfaces()
    const {
  return (m_boundarySurfaces);
}

void Acts::TrackingVolume::connectDenseBoundarySurfaces(
    MutableTrackingVolumeVector& confinedDenseVolumes) {
  if (!confinedDenseVolumes.empty()) {
    NavigationDirection navDir;
    // Walk over each dense volume
    for (auto& confDenseVol : confinedDenseVolumes) {
      // Walk over each boundary surface of the volume
      auto& boundSur = confDenseVol->boundarySurfaces();
      for (unsigned int i = 0; i < boundSur.size(); i++) {
        // Skip empty entries since we do not know the shape of the dense volume
        // and therewith the used indices
        if (boundSur.at(i) == nullptr) {
          continue;
        }

        // Use mother volume as the opposite direction of the already used
        // direction
        auto mutableBs =
            std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(
                boundSur.at(i));
        if (mutableBs->m_oppositeVolume != nullptr &&
            mutableBs->m_alongVolume == nullptr) {
          navDir = forward;
          mutableBs->attachVolume(this, navDir);
        } else {
          if (mutableBs->m_oppositeVolume == nullptr &&
              mutableBs->m_alongVolume != nullptr) {
            navDir = backward;
            mutableBs->attachVolume(this, navDir);
          }
        }

        // Update the boundary
        confDenseVol->updateBoundarySurface((BoundarySurfaceFace)i, mutableBs);
      }
      // Store the volume
      m_confinedDenseVolumes.push_back(std::move(confDenseVol));
    }
  }
}

void Acts::TrackingVolume::createBoundarySurfaces() {
  using Boundary = BoundarySurfaceT<TrackingVolume>;

  // Transform Surfaces To BoundarySurfaces
  auto orientedSurfaces =
      Volume::volumeBounds().orientedSurfaces(m_transform.get());

  m_boundarySurfaces.reserve(orientedSurfaces.size());
  for (auto& osf : orientedSurfaces) {
    TrackingVolume* opposite = nullptr;
    TrackingVolume* along = nullptr;
    if (osf.second == backward) {
      opposite = this;
    } else {
      along = this;
    }
    m_boundarySurfaces.push_back(std::make_shared<const Boundary>(
        std::move(osf.first), opposite, along));
  }
}

void Acts::TrackingVolume::glueTrackingVolume(const GeometryContext& gctx,
                                              BoundarySurfaceFace bsfMine,
                                              TrackingVolume* neighbor,
                                              BoundarySurfaceFace bsfNeighbor) {
  // Find the connection of the two tracking volumes: binR returns the center
  // except for cylindrical volumes
  Vector3D bPosition(binningPosition(gctx, binR));
  Vector3D distance =
      Vector3D(neighbor->binningPosition(gctx, binR) - bPosition);
  // glue to the face
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine =
      boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3D nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  NavigationDirection navDir =
      (nvector.dot(distance) > 0.) ? forward : backward;
  // The easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr) ||
      m_glueVolumeDescriptor->glueVolumes(bsfMine) == nullptr) {
    // the boundary orientation
    auto mutableBSurfaceMine =
        std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(bSurfaceMine);
    mutableBSurfaceMine->attachVolume(neighbor, navDir);
    // Make sure you keep the boundary material if there
    const Surface& neighborSurface =
        neighbor->m_boundarySurfaces.at(bsfNeighbor)->surfaceRepresentation();
    auto neighborMaterial = neighborSurface.surfaceMaterialSharedPtr();
    const Surface& mySurface = bSurfaceMine->surfaceRepresentation();
    auto myMaterial = mySurface.surfaceMaterialSharedPtr();
    // Keep the neighbor material
    if (myMaterial == nullptr and neighborMaterial != nullptr) {
      Surface* myMutbableSurface = const_cast<Surface*>(&mySurface);
      myMutbableSurface->assignSurfaceMaterial(neighborMaterial);
    }
    // Now set it to the neighbor volume
    (neighbor->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
  }
}

void Acts::TrackingVolume::glueTrackingVolumes(
    const GeometryContext& gctx, BoundarySurfaceFace bsfMine,
    const std::shared_ptr<TrackingVolumeArray>& neighbors,
    BoundarySurfaceFace bsfNeighbor) {
  // find the connection of the two tracking volumes : binR returns the center
  // except for cylindrical volumes
  std::shared_ptr<const TrackingVolume> nRefVolume =
      neighbors->arrayObjects().at(0);
  // get the distance
  Vector3D bPosition(binningPosition(gctx, binR));
  Vector3D distance =
      Vector3D(nRefVolume->binningPosition(gctx, binR) - bPosition);
  // take the normal at the binning positio
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine =
      boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3D nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  NavigationDirection navDir =
      (nvector.dot(distance) > 0.) ? forward : backward;
  // the easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr) ||
      !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    auto mutableBSurfaceMine =
        std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(bSurfaceMine);
    mutableBSurfaceMine->attachVolumeArray(neighbors, navDir);
    // now set it to the neighbor volumes - the optised way
    for (auto& nVolume : neighbors->arrayObjects()) {
      auto mutableNVolume = std::const_pointer_cast<TrackingVolume>(nVolume);
      (mutableNVolume->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
    }
  }
}

void Acts::TrackingVolume::assignBoundaryMaterial(
    std::shared_ptr<const ISurfaceMaterial> surfaceMaterial,
    BoundarySurfaceFace bsFace) {
  auto bSurface = m_boundarySurfaces.at(bsFace);
  Surface* surface = const_cast<Surface*>(&bSurface->surfaceRepresentation());
  surface->assignSurfaceMaterial(std::move(surfaceMaterial));
}

void Acts::TrackingVolume::updateBoundarySurface(
    BoundarySurfaceFace bsf,
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs,
    bool checkmaterial) {
  if (checkmaterial) {
    auto cMaterialPtr = m_boundarySurfaces.at(bsf)
                            ->surfaceRepresentation()
                            .surfaceMaterialSharedPtr();
    auto bsMaterial = bs->surfaceRepresentation().surfaceMaterial();
    if (cMaterialPtr != nullptr && bsMaterial == nullptr) {
      Surface* surface = const_cast<Surface*>(&bs->surfaceRepresentation());
      surface->assignSurfaceMaterial(cMaterialPtr);
    }
  }
  m_boundarySurfaces.at(bsf) = std::move(bs);
}

void Acts::TrackingVolume::registerGlueVolumeDescriptor(
    GlueVolumesDescriptor* gvd) {
  delete m_glueVolumeDescriptor;
  m_glueVolumeDescriptor = gvd;
}

Acts::GlueVolumesDescriptor& Acts::TrackingVolume::glueVolumesDescriptor() {
  if (m_glueVolumeDescriptor == nullptr) {
    m_glueVolumeDescriptor = new GlueVolumesDescriptor;
  }
  return (*m_glueVolumeDescriptor);
}

void Acts::TrackingVolume::synchronizeLayers(double envelope) const {
  // case a : Layers exist
  // msgstream << MSG::VERBOSE << "  -> synchronizing Layer dimensions of
  // TrackingVolume '" << volumeName() << "'." << endreq;

  if (m_confinedLayers) {
    // msgstream << MSG::VERBOSE << "  ---> working on " <<
    // m_confinedLayers->arrayObjects().size() << " (material+navigation)
    // layers." << endreq;
    for (auto& clayIter : m_confinedLayers->arrayObjects()) {
      if (clayIter) {
        // @todo implement syncrhonize layer
        //  if (clayIter->surfaceRepresentation().type() == Surface::Cylinder &&
        //  !(center().isApprox(clayIter->surfaceRepresentation().center())) )
        //      clayIter->resizeAndRepositionLayer(volumeBounds(),center(),envelope);
        //  else
        //      clayIter->resizeLayer(volumeBounds(),envelope);
      }  // else
      // msgstream << MSG::WARNING << "  ---> found 0 pointer to layer,
      // indicates problem." << endreq;
    }
  }

  // case b : container volume -> step down
  if (m_confinedVolumes) {
    // msgstream << MSG::VERBOSE << "  ---> no confined layers, working on " <<
    // m_confinedVolumes->arrayObjects().size() << " confined volumes." <<
    // endreq;
    for (auto& cVolumesIter : m_confinedVolumes->arrayObjects()) {
      cVolumesIter->synchronizeLayers(envelope);
    }
  }
}

void Acts::TrackingVolume::interlinkLayers() {
  if (m_confinedLayers) {
    auto& layers = m_confinedLayers->arrayObjects();

    // forward register the last one as the previous one
    //  first <- | -> second, first <- | -> second, first <- | -> second
    const Layer* lastLayer = nullptr;
    for (auto& layerPtr : layers) {
      // we'll need to mutate our confined layers to perform this operation
      Layer& mutableLayer = *(std::const_pointer_cast<Layer>(layerPtr));
      // register the layers
      mutableLayer.m_nextLayerUtility = m_confinedLayers->binUtility();
      mutableLayer.m_nextLayers.first = lastLayer;
      // register the volume
      mutableLayer.encloseTrackingVolume(*this);
      // remember the last layer
      lastLayer = &mutableLayer;
    }
    // backward loop
    lastLayer = nullptr;
    for (auto layerIter = layers.rbegin(); layerIter != layers.rend();
         ++layerIter) {
      // set the other next volume
      Layer& mutableLayer = *(std::const_pointer_cast<Layer>(*layerIter));
      mutableLayer.m_nextLayers.second = lastLayer;
      lastLayer = &mutableLayer;
    }
  }
}

void Acts::TrackingVolume::closeGeometry(
    const IMaterialDecorator* materialDecorator,
    std::map<std::string, const TrackingVolume*>& volumeMap, size_t& vol) {
  // insert the volume into the map
  volumeMap[volumeName()] = this;

  // we can construct the volume ID from this
  auto volumeID = GeometryID().setVolume(++vol);
  // assign the Volume ID to the volume itself
  auto thisVolume = const_cast<TrackingVolume*>(this);
  thisVolume->assignGeoID(volumeID);

  // assign the material if you have a decorator
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*thisVolume);
  }
  if (thisVolume->volumeMaterial() == nullptr && thisVolume->motherVolume() &&
      thisVolume->motherVolume()->volumeMaterial() != nullptr) {
    auto protoMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(
        thisVolume->motherVolume()->volumeMaterial());
    if (protoMaterial == nullptr) {
      thisVolume->assignVolumeMaterial(
          thisVolume->motherVolume()->volumeMaterialSharedPtr());
    }
  }

  this->assignGeoID(volumeID);
  // loop over the boundary surfaces
  GeometryID::Value iboundary = 0;
  // loop over the boundary surfaces
  for (auto& bSurfIter : boundarySurfaces()) {
    // get the intersection soltuion
    auto& bSurface = bSurfIter->surfaceRepresentation();
    // create the boundary surface id
    auto boundaryID = GeometryID(volumeID).setBoundary(++iboundary);
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<Surface*>(&bSurface));
    mutableBSurface.assignGeoID(boundaryID);
    // assign the material if you have a decorator
    if (materialDecorator != nullptr) {
      materialDecorator->decorate(mutableBSurface);
    }
  }

  // A) this is NOT a container volume, volumeID is already incremented
  if (!m_confinedVolumes) {
    // loop over the confined layers
    if (m_confinedLayers) {
      GeometryID::Value ilayer = 0;
      // loop over the layers
      for (auto& layerPtr : m_confinedLayers->arrayObjects()) {
        // create the layer identification
        auto layerID = GeometryID(volumeID).setLayer(++ilayer);
        // now close the geometry
        auto mutableLayerPtr = std::const_pointer_cast<Layer>(layerPtr);
        mutableLayerPtr->closeGeometry(materialDecorator, layerID);
      }
    } else if (m_bvhTop != nullptr) {
      GeometryID::Value isurface = 0;
      for (const auto& descVol : m_descendantVolumes) {
        // Attempt to cast to AbstractVolume: only one we'll handle
        const AbstractVolume* avol =
            dynamic_cast<const AbstractVolume*>(descVol.get());
        if (avol != nullptr) {
          const auto& bndSrf = avol->boundarySurfaces();
          for (const auto& bnd : bndSrf) {
            const auto& srf = bnd->surfaceRepresentation();
            Surface* mutableSurfcePtr = const_cast<Surface*>(&srf);
            auto geoID = GeometryID(volumeID).setSensitive(++isurface);
            mutableSurfcePtr->assignGeoID(geoID);
          }
        }
      }
    }
  } else {
    // B) this is a container volume, go through sub volume
    // do the loop
    for (auto& volumesIter : m_confinedVolumes->arrayObjects()) {
      auto mutableVolumesIter =
          std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->setMotherVolume(this);
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol);
    }
  }

  if (!m_confinedDenseVolumes.empty()) {
    for (auto& volumesIter : m_confinedDenseVolumes) {
      auto mutableVolumesIter =
          std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->setMotherVolume(this);
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol);
    }
  }
}

void Acts::TrackingVolume::visitSurfaces(
    const std::function<void(const Acts::Surface*)>& visitor) const {
  if (!m_confinedVolumes) {
    // no sub volumes => loop over the confined layers
    if (m_confinedLayers) {
      for (const auto& layer : m_confinedLayers->arrayObjects()) {
        if (layer->surfaceArray() == nullptr) {
          // no surface array (?)
          continue;
        }
        for (const auto& srf : layer->surfaceArray()->surfaces()) {
          visitor(srf);
        }
      }
    }
  } else {
    // contains sub volumes
    for (const auto& volume : m_confinedVolumes->arrayObjects()) {
      volume->visitSurfaces(visitor);
    }
  }
}
