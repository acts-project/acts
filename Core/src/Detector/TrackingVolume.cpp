// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include <functional>
#include <utility>

#include "Acts/Detector/DetachedTrackingVolume.hpp"
#include "Acts/Detector/GlueVolumesDescriptor.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

Acts::TrackingVolume::TrackingVolume()
  : Volume()
  , m_material(std::make_shared<const Material>())
  , m_boundarySurfaces()
  , m_confinedLayers(nullptr)
  , m_confinedVolumes(nullptr)
  , m_confinedDetachedVolumes()
  , m_confinedDenseVolumes()
  , m_confinedArbitraryLayers()
  , m_name("undefined")
{
}

Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<const Transform3D>                htrans,
    VolumeBoundsPtr                                   volbounds,
    const std::shared_ptr<const TrackingVolumeArray>& containedVolumeArray,
    const std::string&                                volumeName)
  : Volume(std::move(htrans), std::move(volbounds))
  , m_material(std::make_shared<const Material>())
  , m_boundarySurfaces()
  , m_confinedLayers(nullptr)
  , m_confinedVolumes(containedVolumeArray)
  , m_confinedDetachedVolumes()
  , m_confinedDenseVolumes()
  , m_confinedArbitraryLayers()
  , m_name(volumeName)
{
  createBoundarySurfaces();
  interlinkLayers();
}

// constructor for arguments
Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<const Transform3D>         htrans,
    VolumeBoundsPtr                            volbounds,
    std::shared_ptr<const Material>            matprop,
    std::unique_ptr<const LayerArray>          staticLayerArray,
    const LayerVector&                         arbitraryLayerVector,
    std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
    const TrackingVolumeVector&                denseVolumeVector,
    const DetachedVolumeVector&                detachedVolumeVector,
    const std::string&                         volumeName)
  : Volume(std::move(htrans), std::move(volbounds))
  , m_material(std::move(matprop))
  , m_confinedLayers(std::move(staticLayerArray))
  , m_confinedVolumes(std::move(containedVolumeArray))
  , m_confinedDetachedVolumes(detachedVolumeVector)
  , m_confinedDenseVolumes(denseVolumeVector)
  , m_confinedArbitraryLayers(arbitraryLayerVector)
  , m_name(volumeName)
{
  createBoundarySurfaces();
  interlinkLayers();
}

Acts::TrackingVolume::~TrackingVolume()
{
  delete m_glueVolumeDescriptor;
}

const Acts::Layer*
Acts::TrackingVolume::associatedLayer(const Vector3D& gp) const
{
  // confined static layers - highest hierarchy
  if (m_confinedLayers) {
    return (m_confinedLayers->object(gp).get());
  }

  // confined arbitrary
  if (!m_confinedArbitraryLayers.empty()) {
    for (auto& layer : m_confinedArbitraryLayers) {
      if (layer->isOnLayer(gp)) {
        return layer.get();
      }
    }
  }

  // return the null pointer
  return nullptr;
}

const Acts::TrackingVolume*
Acts::TrackingVolume::trackingVolume(const Vector3D& gp) const
{
  // confined static volumes - highest hierarchy
  if (m_confinedVolumes) {
    return (m_confinedVolumes->object(gp).get());
  }

  // if no static volumes are there, detached is next hierarchy
  if (!m_confinedDetachedVolumes.empty()) {
    for (auto& detachedVolume : m_confinedDetachedVolumes) {
      if (detachedVolume->trackingVolume()->inside(gp, 0.001)) {
        {
          return detachedVolume->trackingVolume();
        }
      }
    }
  }

  // if no static volumes or detached volumes are there, search for dense
  // volumes
  if (!m_confinedDenseVolumes.empty()) {
    for (auto& denseVolume : m_confinedDenseVolumes) {
      if (denseVolume->inside(gp, 0.001)) {
        return denseVolume.get();
      }
    }
  }

  // there is no lower sub structure
  return this;
}

const Acts::DetachedVolumeVector*
Acts::TrackingVolume::detachedTrackingVolumes(const Vector3D& gp,
                                              double          tol) const
{
  // create a new vector
  DetachedVolumeVector* currVols = new DetachedVolumeVector;
  // get the volumes were the position is inside
  if (!m_confinedDetachedVolumes.empty()) {
    for (auto& detachedVolume : m_confinedDetachedVolumes) {
      if (detachedVolume->trackingVolume()->inside(gp, tol)) {
        currVols->push_back(detachedVolume);
      }
    }
  }
  // return the volumes that are inside
  return currVols;
}

void
Acts::TrackingVolume::sign(GeometrySignature geosign, GeometryType geotype)
{
  // never overwrite what is already signed, that's a crime
  if (m_geometrySignature == Unsigned) {
    m_geometrySignature = geosign;
  }
  m_geometryType = geotype;

  // confined static volumes
  if (m_confinedVolumes) {
    for (auto& volumesIter : (m_confinedVolumes->arrayObjects())) {
      auto mutableVolumesIter
          = std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->sign(geosign, geotype);
    }
  }

  // same procedure for the detached volumes
  if (!m_confinedDetachedVolumes.empty()) {
    for (auto& volumesIter : m_confinedDetachedVolumes) {
      auto mutableVolumesIter
          = std::const_pointer_cast<DetachedTrackingVolume>(volumesIter);
      mutableVolumesIter->sign(geosign, geotype);
    }
  }

  // finally for confined dense volumes
  if (!m_confinedDenseVolumes.empty()) {
    for (auto& volumesIter : m_confinedDenseVolumes) {
      auto mutableVolumesIter
          = std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->sign(geosign, geotype);
    }
  }
}

const std::
    vector<std::shared_ptr<const Acts::BoundarySurfaceT<Acts::TrackingVolume>>>&
    Acts::TrackingVolume::boundarySurfaces() const
{
  return (m_boundarySurfaces);
}

void
Acts::TrackingVolume::createBoundarySurfaces()
{
  // transform Surfaces To BoundarySurfaces
  std::vector<std::shared_ptr<const Surface>> surfaces
      = Volume::volumeBounds().decomposeToSurfaces(m_transform);

  // counter to flip the inner/outer position for Cylinders
  int    sfCounter = 0;
  size_t sfNumber  = surfaces.size();

  for (auto& sf : surfaces) {
    // flip inner/outer for cylinders
    TrackingVolume* inner
        = (sf->type() == Surface::Cylinder && sfCounter == 3 && sfNumber > 3)
        ? nullptr
        : this;
    TrackingVolume* outer = (inner) != nullptr ? nullptr : this;
    // create the boundary surface
    m_boundarySurfaces.push_back(
        std::make_shared<const BoundarySurfaceT<TrackingVolume>>(
            std::move(sf), inner, outer));
    // increase the counter
    ++sfCounter;
  }
}

void
Acts::TrackingVolume::glueTrackingVolume(
    BoundarySurfaceFace                    bsfMine,
    const std::shared_ptr<TrackingVolume>& neighbor,
    BoundarySurfaceFace                    bsfNeighbor)
{
  // find the connection of the two tracking volumes : binR returns the center
  // except for cylindrical volumes
  Vector3D bPosition(binningPosition(binR));
  Vector3D distance = Vector3D(neighbor->binningPosition(binR) - bPosition);
  // glue to the face
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine
      = boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
  // estimate the orientation
  BoundaryOrientation bOrientation
      = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
  // the easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr)
      || !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    auto mutableBSurfaceMine
        = std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(
            bSurfaceMine);
    mutableBSurfaceMine->attachVolume(neighbor, bOrientation);
    // now set it to the neighbor volume - the optised way
    (neighbor->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
  }
}

void
Acts::TrackingVolume::glueTrackingVolumes(
    BoundarySurfaceFace                         bsfMine,
    const std::shared_ptr<TrackingVolumeArray>& neighbors,
    BoundarySurfaceFace                         bsfNeighbor)
{
  // find the connection of the two tracking volumes : binR returns the center
  // except for cylindrical volumes
  std::shared_ptr<const TrackingVolume> nRefVolume
      = neighbors->arrayObjects().at(0);
  // get the distance
  Vector3D bPosition(binningPosition(binR));
  Vector3D distance = Vector3D(nRefVolume->binningPosition(binR) - bPosition);
  // take the normal at the binning positio
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine
      = boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
  // estimate the orientation
  BoundaryOrientation bOrientation
      = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
  // the easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr)
      || !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    auto mutableBSurfaceMine
        = std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(
            bSurfaceMine);
    mutableBSurfaceMine->attachVolumeArray(neighbors, bOrientation);
    // now set it to the neighbor volumes - the optised way
    for (auto& nVolume : neighbors->arrayObjects()) {
      auto mutableNVolume = std::const_pointer_cast<TrackingVolume>(nVolume);
      (mutableNVolume->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
    }
  }
}

void
Acts::TrackingVolume::updateBoundarySurface(
    BoundarySurfaceFace                                     bsf,
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs)
{
  m_boundarySurfaces.at(bsf) = std::move(bs);
}

void
Acts::TrackingVolume::registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd)
{
  delete m_glueVolumeDescriptor;
  m_glueVolumeDescriptor = gvd;
}

Acts::GlueVolumesDescriptor&
Acts::TrackingVolume::glueVolumesDescriptor()
{
  if (m_glueVolumeDescriptor == nullptr) {
    m_glueVolumeDescriptor = new GlueVolumesDescriptor;
  }
  return (*m_glueVolumeDescriptor);
}

void
Acts::TrackingVolume::synchronizeLayers(double envelope) const
{
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

void
Acts::TrackingVolume::interlinkLayers()
{
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
      lastLayer                        = &mutableLayer;
    }
  }
}

void
Acts::TrackingVolume::closeGeometry(
    std::map<std::string, const TrackingVolume*>& volumeMap,
    size_t& vol)
{
  // insert the volume into the map
  volumeMap[volumeName()] = this;

  // we can construct the volume ID from this
  GeometryID volumeID(0);
  volumeID.add(++vol, GeometryID::volume_mask);
  // assign the Volume ID to the volume itself
  auto thisVolume = const_cast<TrackingVolume*>(this);
  thisVolume->assignGeoID(volumeID);

  // loop over the boundary surfaces
  geo_id_value iboundary = 0;
  // loop over the boundary surfaces
  for (auto& bSurfIter : boundarySurfaces()) {
    // get the intersection soltuion
    auto& bSurface = bSurfIter->surfaceRepresentation();
    // create the boundary surface id
    GeometryID boundaryID = volumeID;
    boundaryID.add(++iboundary, GeometryID::boundary_mask);
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<Surface*>(&bSurface));
    mutableBSurface.assignGeoID(boundaryID);
  }

  // A) this is NOT a container volume, volumeID is already incremented
  if (!m_confinedVolumes) {
    // loop over the confined layers
    if (m_confinedLayers) {
      geo_id_value ilayer = 0;
      // loop over the layers
      for (auto& layerPtr : m_confinedLayers->arrayObjects()) {
        // create the layer identification
        GeometryID layerID = volumeID;
        layerID.add(++ilayer, GeometryID::layer_mask);
        // now close the geometry
        auto mutableLayerPtr = std::const_pointer_cast<Layer>(layerPtr);
        mutableLayerPtr->closeGeometry(layerID);
      }
    }
  } else {
    // B) this is a container volume, go through sub volume
    // do the loop
    for (auto& volumesIter : m_confinedVolumes->arrayObjects()) {
      auto mutableVolumesIter
          = std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->closeGeometry(volumeMap, vol);
    }
  }

  // @todo update that
  // auto confinedDenseVolumes= tvol.confinedDenseVolumes();
  // if (!confinedDenseVolumes.empty()) {
  //   for (auto& volumesIter : confinedDenseVolumes)
  //     if (volumesIter) closeGeometry(*volumesIter, &tvol, ++cCounter);
  // }
  //
  // // should detached tracking volumes be part of the tracking geometry ? */
  // auto confinedDetachedVolumes = tvol.confinedDetachedVolumes();
  // if (!confinedDetachedVolumes.empty()) {
  //   for (auto& volumesIter : confinedDetachedVolumes)
  //     if (volumesIter
  //         && tvol.inside(volumesIter->trackingVolume()->center(), 0.))
  //       closeGeometry(*(volumesIter->trackingVolume()), &tvol, ++cCounter);
  // }
  //
}

void
Acts::TrackingVolume::visitSurfaces(
    const std::function<void(const Acts::Surface*)>& visitor) const
{
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
