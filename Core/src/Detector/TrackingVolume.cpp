// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Detector/DetachedTrackingVolume.hpp"
#include "ACTS/Detector/GlueVolumesDescriptor.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"

Acts::TrackingVolume::TrackingVolume()
  : Volume()
  , m_material(std::make_shared<Acts::Material>())
  , m_motherVolume()
  , m_boundarySurfaces()
  , m_confinedLayers(nullptr)
  , m_confinedVolumes(nullptr)
  , m_confinedDetachedVolumes()
  , m_confinedDenseVolumes()
  , m_confinedArbitraryLayers()
  , m_glueVolumeDescriptor(nullptr)
  , m_geometrySignature(Unsigned)
  , m_geometryType(NumberOfGeometryTypes)
  , m_name("undefined")
  , m_colorCode(20)
{
}

Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<Transform3D>                     htrans,
    VolumeBoundsPtr                                  volbounds,
    const std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
    const std::string&                               volumeName)
  : Volume(htrans, volbounds)
  , m_material(std::make_shared<Acts::Material>())
  , m_motherVolume(nullptr)
  , m_boundarySurfaces()
  , m_confinedLayers(nullptr)
  , m_confinedVolumes(containedVolumeArray)
  , m_confinedDetachedVolumes()
  , m_confinedDenseVolumes()
  , m_confinedArbitraryLayers()
  , m_glueVolumeDescriptor(nullptr)
  , m_geometrySignature(Unsigned)
  , m_geometryType(NumberOfGeometryTypes)
  , m_name(volumeName)
  , m_colorCode(20)
{
  createBoundarySurfaces();
  interlinkLayers();
}

// constructor for arguments
Acts::TrackingVolume::TrackingVolume(
    std::shared_ptr<Transform3D>               htrans,
    VolumeBoundsPtr                            volbounds,
    std::shared_ptr<Material>                  matprop,
    std::unique_ptr<const LayerArray>          staticLayerArray,
    const LayerVector                          arbitraryLayerVector,
    std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
    const TrackingVolumeVector                 denseVolumeVector,
    const DetachedVolumeVector                 detachedVolumeVector,
    const std::string&                         volumeName)
  : Volume(htrans, volbounds)
  , m_material(matprop)
  , m_motherVolume(nullptr)
  , m_confinedLayers(std::move(staticLayerArray))
  , m_confinedVolumes(containedVolumeArray)
  , m_confinedDetachedVolumes(detachedVolumeVector)
  , m_confinedDenseVolumes(denseVolumeVector)
  , m_confinedArbitraryLayers(arbitraryLayerVector)
  , m_glueVolumeDescriptor(nullptr)
  , m_geometrySignature(Unsigned)
  , m_geometryType(NumberOfGeometryTypes)
  , m_name(volumeName)
  , m_colorCode(20)
{
  createBoundarySurfaces();
  interlinkLayers();
}

Acts::TrackingVolume::TrackingVolume(const TrackingVolume& tvol,
                                     const Transform3D&    shift,
                                     const std::string&    volumeName)
  : Volume(tvol, &shift)
  , m_material(tvol.m_material)
  , m_motherVolume(tvol.motherVolume())
  , m_confinedLayers(nullptr)
  , m_confinedVolumes(nullptr)
  , m_confinedDetachedVolumes()
  , m_confinedDenseVolumes()
  , m_confinedArbitraryLayers()
  , m_glueVolumeDescriptor(nullptr)
  , m_geometrySignature(tvol.geometrySignature())
  , m_geometryType(tvol.geometryType())
  , m_name(volumeName)
  , m_colorCode(20)
{
  //< @TODO implement - requires cloneWithShift for BinUtility and an
  // orderPosition() addon to GeometryObjects
}

Acts::TrackingVolume::~TrackingVolume()
{
  delete m_glueVolumeDescriptor;
}

const Acts::Layer*
Acts::TrackingVolume::associatedLayer(const Vector3D& gp) const
{
  // confined static layers - highest hierarchy
  if (m_confinedLayers) return (m_confinedLayers->object(gp).get());

  // confined arbitrary
  if (!m_confinedArbitraryLayers.empty())
    for (auto& layer : m_confinedArbitraryLayers)
      if (layer->isOnLayer(gp)) return layer.get();

  // return the null pointer
  return nullptr;
}

const Acts::TrackingVolume*
Acts::TrackingVolume::trackingVolume(const Vector3D& gp) const
{
  // confined static volumes - highest hierarchy
  if (m_confinedVolumes) return (m_confinedVolumes->object(gp).get());

  // if no static volumes are there, detached is next hierarchy
  if (!m_confinedDetachedVolumes.empty())
    for (auto& detachedVolume : m_confinedDetachedVolumes)
      if (detachedVolume->trackingVolume()->inside(gp, 0.001))
        return detachedVolume->trackingVolume();

  // if no static volumes or detached volumes are there, search for dense
  // volumes
  if (!m_confinedDenseVolumes.empty())
    for (auto& denseVolume : m_confinedDenseVolumes)
      if (denseVolume->inside(gp, 0.001)) return denseVolume.get();

  // there is no lower sub structure
  return this;
}

const Acts::TrackingVolume*
Acts::TrackingVolume::nextVolume(const Vector3D& gp,
                                 const Vector3D& dir,
                                 PropDirection   pDir) const
{
  // get the boundary surfaces & intersect them
  const TrackingVolume* nVolume = 0;
  // fix the direction once
  bool     forceDir   = (pDir == alongMomentum || pDir == oppositeMomentum);
  double   dirScalor  = (pDir == oppositeMomentum) ? -1. : 1.;
  Vector3D cDir       = dirScalor * dir;
  double   pathLength = 10e10;
  // now loop through the and find the closest
  auto bSurfaces = boundarySurfaces();
  for (auto& bSurfIter : bSurfaces) {
    // get the intersection soltuion
    Intersection sfI = bSurfIter->surfaceRepresentation().intersectionEstimate(
        gp, cDir, forceDir, true);
    if (sfI.valid
        && (sfI.pathLength * sfI.pathLength) < (pathLength * pathLength)) {
      // assign the next Volume
      PropDirection attachedDir
          = sfI.pathLength > 0. ? alongMomentum : oppositeMomentum;
      pathLength = sfI.pathLength;
      nVolume    = bSurfIter->attachedVolume(gp, cDir, attachedDir);
    }
  }
  return nVolume;
}

const Acts::DetachedVolumeVector*
Acts::TrackingVolume::detachedTrackingVolumes(const Vector3D& gp,
                                              double          tol) const
{
  // create a new vector
  DetachedVolumeVector* currVols = new DetachedVolumeVector;
  // get the volumes were the position is inside
  if (!m_confinedDetachedVolumes.empty())
    for (auto& detachedVolume : m_confinedDetachedVolumes)
      if (detachedVolume->trackingVolume()->inside(gp, tol))
        currVols->push_back(detachedVolume);
  // return the volumes that are inside
  return currVols;
}

void
Acts::TrackingVolume::addMaterial(std::shared_ptr<const Material> mprop,
                                  float                           fact)
{
  // assume the scaling factor refers to the volume scaling
  float flin = pow(fact, 0.33);
  // average X0
  double invX0     = m_material->X0 > 0. ? 1. / m_material->X0 : 0.;
  double sum_invX0 = invX0 + flin / mprop->X0;
  float  X0        = 1. / sum_invX0;
  // average L0
  double invL0     = m_material->L0 > 0. ? 1. / m_material->L0 : 0.;
  double sum_invL0 = invL0 + flin / mprop->L0;
  float  L0        = 1. / sum_invL0;
  // add density
  float rho1 = m_material->rho;
  float rho  = rho1 + fact * mprop->rho;
  // averageZ
  float n1 = m_material->Z > 0. ? rho1 / m_material->Z : 0.;
  float n2 = fact * mprop->rho / mprop->Z;
  float Z  = rho / (n1 + n2);
  // averageA
  n1      = m_material->A > 0. ? rho1 / m_material->A : 0.;
  n2      = fact * mprop->rho / mprop->A;
  float A = rho / (n1 + n2);
  // mean energy loss (linear scaling)
  float dEdX = m_material->dEdX + flin * mprop->dEdX;

  m_material.reset(new Material(X0, L0, A, Z, rho, dEdX));
  // m_material.reset(std::make_shared<Material>(X0,L0,A,Z,rho,dEdX));
}

void
Acts::TrackingVolume::sign(GeometrySignature geosign,
                           GeometryType      geotype) const
{
  // never overwrite what is already signed, that's a crime
  if (m_geometrySignature == Unsigned) m_geometrySignature = geosign;
  m_geometryType                                           = geotype;

  // confined static volumes
  if (m_confinedVolumes)
    for (auto& volumesIter : (m_confinedVolumes->arrayObjects()))
      volumesIter->sign(geosign, geotype);

  // same procedure for the detached volumes
  if (!m_confinedDetachedVolumes.empty())
    for (auto& volumesIter : m_confinedDetachedVolumes)
      volumesIter->sign(geosign, geotype);

  // finally for confined dense volumes
  if (!m_confinedDenseVolumes.empty())
    for (auto& volumesIter : m_confinedDenseVolumes)
      volumesIter->sign(geosign, geotype);
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
  const std::vector<const Surface*> surfaces
      = Volume::volumeBounds().decomposeToSurfaces(m_transform);

  // counter to flip the inner/outer position for Cylinders
  int sfCounter = 0;
  size_t sfNumber  = surfaces.size();

  for (auto& sf : surfaces) {
    // flip inner/outer for cylinders
    TrackingVolume* inner
        = (sf->type() == Surface::Cylinder && sfCounter == 3 && sfNumber > 3)
        ? nullptr
        : this;
    TrackingVolume* outer = (inner) ? nullptr : this;
    // create the boundary surface
    m_boundarySurfaces.push_back(
        std::make_shared<const BoundarySurfaceT<TrackingVolume>>(
            std::unique_ptr<const Surface>(sf), inner, outer));
    // increase the counter
    ++sfCounter;
  }
}

void
Acts::TrackingVolume::glueTrackingVolume(
    BoundarySurfaceFace                   bsfMine,
    std::shared_ptr<const TrackingVolume> neighbor,
    BoundarySurfaceFace                   bsfNeighbor) const
{
  // find the connection of the two tracking volumes : binR returns the center
  // except for cylindrical volumes
  Vector3D bPosition(binningPosition(binR));
  Vector3D distance = Vector3D(neighbor->binningPosition(binR) - bPosition);
  // glue to the face
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine
      = boundarySurfaces().at(bsfMine);
  // @TODO - complex glueing could be possible
  // with actual intersection for the normal vector
  Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
  // estimate the orientation
  BoundaryOrientation bOrientation
      = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
  // the easy case :
  // - no glue volume descriptors on either side
  if (!m_glueVolumeDescriptor
      || !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    bSurfaceMine->attachVolume(neighbor, bOrientation);
    // now set it to the neighbor volume - the optised way
    (neighbor->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
  }
}

void
Acts::TrackingVolume::glueTrackingVolumes(
    BoundarySurfaceFace                        bsfMine,
    std::shared_ptr<const TrackingVolumeArray> neighbors,
    BoundarySurfaceFace                        bsfNeighbor) const
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
  // @TODO - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3D normal = bSurfaceMine->surfaceRepresentation().normal(bPosition);
  // estimate the orientation
  BoundaryOrientation bOrientation
      = (normal.dot(distance) > 0.) ? outsideVolume : insideVolume;
  // the easy case :
  // - no glue volume descriptors on either side
  if (!m_glueVolumeDescriptor
      || !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    bSurfaceMine->attachVolumeArray(neighbors, bOrientation);
    // now set it to the neighbor volumes - the optised way
    for (auto& nVolume : neighbors->arrayObjects())
      (nVolume->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
  }
}

void
Acts::TrackingVolume::updateBoundarySurface(
    BoundarySurfaceFace                                     bsf,
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs) const
{
  m_boundarySurfaces.at(bsf) = bs;
}

void
Acts::TrackingVolume::registerGlueVolumeDescriptor(
    GlueVolumesDescriptor* gvd) const
{
  delete m_glueVolumeDescriptor;
  m_glueVolumeDescriptor = gvd;
}

const Acts::GlueVolumesDescriptor&
Acts::TrackingVolume::glueVolumesDescriptor() const
{
  if (!m_glueVolumeDescriptor)
    m_glueVolumeDescriptor = new GlueVolumesDescriptor;
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
    for (auto& clayIter : m_confinedLayers->arrayObjects())
      if (clayIter) {
        // @TODO implement syncrhonize layer
        //  if (clayIter->surfaceRepresentation().type() == Surface::Cylinder &&
        //  !(center().isApprox(clayIter->surfaceRepresentation().center())) )
        //      clayIter->resizeAndRepositionLayer(volumeBounds(),center(),envelope);
        //  else
        //      clayIter->resizeLayer(volumeBounds(),envelope);
      }  // else
    // msgstream << MSG::WARNING << "  ---> found 0 pointer to layer, indicates
    // problem." << endreq;
  }

  // case b : container volume -> step down
  if (m_confinedVolumes) {
    // msgstream << MSG::VERBOSE << "  ---> no confined layers, working on " <<
    // m_confinedVolumes->arrayObjects().size() << " confined volumes." <<
    // endreq;
    for (auto& cVolumesIter : m_confinedVolumes->arrayObjects())
      cVolumesIter->synchronizeLayers(envelope);
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
      // register the layers
      (*layerPtr).m_nextLayerUtility = m_confinedLayers->binUtility();
      (*layerPtr).m_nextLayers.first = lastLayer;
      // register the volume
      (*layerPtr).encloseTrackingVolume(*this);
      // remember the last layer
      lastLayer = layerPtr.get();
    }
    // backward loop
    lastLayer = nullptr;
    for (auto layerIter = layers.rbegin(); layerIter != layers.rend();
         ++layerIter) {
      // set the other next volume
      (**layerIter).m_nextLayers.second = lastLayer;
      lastLayer                         = (*layerIter).get();
    }
  }
}

void
Acts::TrackingVolume::closeGeometry(
    const GeometryID& volumeID,
    std::map<std::string, const TrackingVolume*>& volumeMap) const
{
  // insert the volume into the map
  volumeMap[volumeName()] = this;

  // A) this is NOT a container volume, volumeID is already incremented
  if (!m_confinedVolumes) {
    // assign the Volume ID to the volume itself
    assignGeoID(volumeID);
    // loop over the boundary surfaces
    geo_id_value iboundary = 0;
    for (auto& bSurfIter : boundarySurfaces()) {
      // get the intersection soltuion
      auto& bSurface = bSurfIter->surfaceRepresentation();
      // create the boundary surface id
      GeometryID boundaryID = volumeID;
      boundaryID += (++iboundary << GeometryID::boundary_shift);
      // now assign to the boundary surface
      bSurface.assignGeoID(boundaryID);
    }

    // loop over the confined layers
    if (m_confinedLayers) {
      geo_id_value ilayer = 0;
      // loop over the layers
      for (auto& layerPtr : m_confinedLayers->arrayObjects()) {
        // create the layer identification
        GeometryID layerID = volumeID;
        layerID += (++ilayer << GeometryID::layer_shift);
        // now close the geometry
        layerPtr->closeGeometry(layerID);
      }
    }
  } else {
    // B) this is a container volume, go through sub volume
    // the counter upwards
    geo_id_value ivolume = 0;
    // do the loop
    for (auto& volumesIter : m_confinedVolumes->arrayObjects()) {
      GeometryID currentID = volumeID;
      // only increase the counter if it's not a container volume
      if (!volumesIter->confinedVolumes()) {
        /// we count the volume ID up
        currentID += (++ivolume << GeometryID::volume_shift);
      }
      volumesIter->closeGeometry(currentID, volumeMap);
    }
  }

  // @TODO update that
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
