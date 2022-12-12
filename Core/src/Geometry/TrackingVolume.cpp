// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingVolume.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GlueVolumesDescriptor.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Ray.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <string>
#include <utility>

Acts::TrackingVolume::TrackingVolume(
    const Transform3& transform, VolumeBoundsPtr volbounds,
    const std::shared_ptr<const TrackingVolumeArray>& containedVolumeArray,
    const std::string& volumeName)
    : Volume(transform, std::move(volbounds)),
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
    const Transform3& transform, VolumeBoundsPtr volumeBounds,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    std::unique_ptr<const LayerArray> staticLayerArray,
    std::shared_ptr<const TrackingVolumeArray> containedVolumeArray,
    MutableTrackingVolumeVector denseVolumeVector,
    const std::string& volumeName)
    : Volume(transform, std::move(volumeBounds)),
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
    const Transform3& transform, VolumeBoundsPtr volbounds,
    std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
    std::vector<std::unique_ptr<const Volume>> descendants,
    const Volume::BoundingBox* top,
    std::shared_ptr<const IVolumeMaterial> volumeMaterial,
    const std::string& volumeName)
    : Volume(transform, std::move(volbounds)),
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
    const GeometryContext& /*gctx*/, const Vector3& position,
    const double tol) const {
  // confined static volumes - highest hierarchy
  if (m_confinedVolumes) {
    return (m_confinedVolumes->object(position).get());
  }

  // search for dense volumes
  if (!m_confinedDenseVolumes.empty()) {
    for (auto& denseVolume : m_confinedDenseVolumes) {
      if (denseVolume->inside(position, tol)) {
        return denseVolume.get();
      }
    }
  }

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
    NavigationDirection navDir = NavigationDirection::Forward;
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
          navDir = NavigationDirection::Forward;
          mutableBs->attachVolume(this, navDir);
        } else {
          if (mutableBs->m_oppositeVolume == nullptr &&
              mutableBs->m_alongVolume != nullptr) {
            navDir = NavigationDirection::Backward;
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
  auto orientedSurfaces = Volume::volumeBounds().orientedSurfaces(m_transform);

  m_boundarySurfaces.reserve(orientedSurfaces.size());
  for (auto& osf : orientedSurfaces) {
    TrackingVolume* opposite = nullptr;
    TrackingVolume* along = nullptr;
    if (osf.second == NavigationDirection::Backward) {
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
  Vector3 bPosition(binningPosition(gctx, binR));
  Vector3 distance = Vector3(neighbor->binningPosition(gctx, binR) - bPosition);
  // glue to the face
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine =
      boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3 nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  NavigationDirection navDir = (nvector.dot(distance) > 0.)
                                   ? NavigationDirection::Forward
                                   : NavigationDirection::Backward;
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
  Vector3 bPosition(binningPosition(gctx, binR));
  Vector3 distance =
      Vector3(nRefVolume->binningPosition(gctx, binR) - bPosition);
  // take the normal at the binning positio
  std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bSurfaceMine =
      boundarySurfaces().at(bsfMine);
  // @todo - complex glueing could be possible with actual intersection for the
  // normal vector
  Vector3 nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  NavigationDirection navDir = (nvector.dot(distance) > 0.)
                                   ? NavigationDirection::Forward
                                   : NavigationDirection::Backward;
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
    std::unordered_map<GeometryIdentifier, const TrackingVolume*>& volumeMap,
    size_t& vol, const GeometryIdentifierHook& hook) {
  // we can construct the volume ID from this
  auto volumeID = GeometryIdentifier().setVolume(++vol);
  // assign the Volume ID to the volume itself
  auto thisVolume = const_cast<TrackingVolume*>(this);
  thisVolume->assignGeometryId(volumeID);

  // insert the volume into the map
  volumeMap[volumeID] = thisVolume;

  // assign the material if you have a decorator
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*thisVolume);
  }
  if (thisVolume->volumeMaterial() == nullptr &&
      thisVolume->motherVolume() != nullptr &&
      thisVolume->motherVolume()->volumeMaterial() != nullptr) {
    auto protoMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(
        thisVolume->motherVolume()->volumeMaterial());
    if (protoMaterial == nullptr) {
      thisVolume->assignVolumeMaterial(
          thisVolume->motherVolume()->volumeMaterialSharedPtr());
    }
  }

  this->assignGeometryId(volumeID);
  // loop over the boundary surfaces
  GeometryIdentifier::Value iboundary = 0;
  // loop over the boundary surfaces
  for (auto& bSurfIter : boundarySurfaces()) {
    // get the intersection soltuion
    auto& bSurface = bSurfIter->surfaceRepresentation();
    // create the boundary surface id
    auto boundaryID = GeometryIdentifier(volumeID).setBoundary(++iboundary);
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<Surface*>(&bSurface));
    mutableBSurface.assignGeometryId(boundaryID);
    // Assigne material if you have a decorator
    if (materialDecorator != nullptr) {
      materialDecorator->decorate(mutableBSurface);
    }
  }

  // A) this is NOT a container volume, volumeID is already incremented
  if (!m_confinedVolumes) {
    // loop over the confined layers
    if (m_confinedLayers) {
      GeometryIdentifier::Value ilayer = 0;
      // loop over the layers
      for (auto& layerPtr : m_confinedLayers->arrayObjects()) {
        // create the layer identification
        auto layerID = GeometryIdentifier(volumeID).setLayer(++ilayer);
        // now close the geometry
        auto mutableLayerPtr = std::const_pointer_cast<Layer>(layerPtr);
        mutableLayerPtr->closeGeometry(materialDecorator, layerID, hook);
      }
    } else if (m_bvhTop != nullptr) {
      GeometryIdentifier::Value isurface = 0;
      for (const auto& descVol : m_descendantVolumes) {
        // Attempt to cast to AbstractVolume: only one we'll handle
        const AbstractVolume* avol =
            dynamic_cast<const AbstractVolume*>(descVol.get());
        if (avol != nullptr) {
          const auto& bndSrf = avol->boundarySurfaces();
          for (const auto& bnd : bndSrf) {
            const auto& srf = bnd->surfaceRepresentation();
            Surface* mutableSurfcePtr = const_cast<Surface*>(&srf);
            auto geoID = GeometryIdentifier(volumeID).setSensitive(++isurface);
            mutableSurfcePtr->assignGeometryId(geoID);
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
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol,
                                        hook);
    }
  }

  if (!m_confinedDenseVolumes.empty()) {
    for (auto& volumesIter : m_confinedDenseVolumes) {
      auto mutableVolumesIter =
          std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->setMotherVolume(this);
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol,
                                        hook);
    }
  }
}

// Returns the boundary surfaces ordered in probability to hit them based on
boost::container::small_vector<Acts::BoundaryIntersection, 4>
Acts::TrackingVolume::compatibleBoundaries(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Surface>& options,
    LoggerWrapper logger) const {
  ACTS_VERBOSE("Finding compatibleBoundaries");
  // Loop over boundarySurfaces and calculate the intersection
  auto excludeObject = options.startObject;
  boost::container::small_vector<Acts::BoundaryIntersection, 4> bIntersections;

  // The signed direction: solution (except overstepping) is positive
  auto sDirection = options.navDir * direction;

  // The Limits: current, path & overstepping
  double pLimit = options.pathLimit;
  double oLimit = options.overstepLimit;

  // Helper function to test intersection
  auto checkIntersection =
      [&](SurfaceIntersection& sIntersection,
          const BoundarySurface* bSurface) -> BoundaryIntersection {
    // Avoid doing anything if that's a rotten apple already
    if (!sIntersection) {
      return BoundaryIntersection();
    }

    if (options.forceIntersectBoundaries and
        sIntersection.intersection.pathLength * options.navDir > 0) {
      const bool coCriterion =
          std::abs(sIntersection.intersection.pathLength) < std::abs(oLimit);
      ACTS_VERBOSE("Forcing intersection with surface "
                   << bSurface->surfaceRepresentation().geometryId());
      if (coCriterion) {
        ACTS_VERBOSE("Intersection forced successfully ");
        ACTS_VERBOSE("- intersection path length "
                     << std::abs(sIntersection.intersection.pathLength)
                     << " < overstep limit " << std::abs(oLimit));
        sIntersection.intersection.pathLength *= options.navDir;
        return BoundaryIntersection(sIntersection.intersection, bSurface,
                                    sIntersection.object);
      }
      ACTS_VERBOSE("Can't force intersection: ");
      ACTS_VERBOSE("- intersection path length "
                   << std::abs(sIntersection.intersection.pathLength)
                   << " > overstep limit " << std::abs(oLimit));
    }

    ACTS_VERBOSE("Check intersection with surface "
                 << bSurface->surfaceRepresentation().geometryId());
    if (detail::checkIntersection(sIntersection.intersection, pLimit, oLimit,
                                  s_onSurfaceTolerance, logger)) {
      sIntersection.intersection.pathLength *= options.navDir;
      return BoundaryIntersection(sIntersection.intersection, bSurface,
                                  sIntersection.object);
    }

    if (sIntersection.alternative) {
      ACTS_VERBOSE("Consider alternative");
      if (detail::checkIntersection(sIntersection.alternative, pLimit, oLimit,
                                    s_onSurfaceTolerance, logger)) {
        sIntersection.alternative.pathLength *= options.navDir;
        return BoundaryIntersection(sIntersection.alternative, bSurface,
                                    sIntersection.object);
        ;
      }
    } else {
      ACTS_VERBOSE("No alternative for intersection");
    }

    ACTS_VERBOSE("No intersection accepted");
    return BoundaryIntersection();
  };

  /// Helper function to process boundary surfaces
  auto processBoundaries =
      [&](const TrackingVolumeBoundaries& bSurfaces) -> void {
    ACTS_VERBOSE("Processing boundaries");
    // Loop over the boundary surfaces
    for (auto& bsIter : bSurfaces) {
      // Get the boundary surface pointer
      const auto& bSurfaceRep = bsIter->surfaceRepresentation();
      ACTS_VERBOSE("Consider boundary surface " << bSurfaceRep.geometryId()
                                                << " :\n"
                                                << std::tie(bSurfaceRep, gctx));

      // Exclude the boundary where you are on
      if (excludeObject != &bSurfaceRep) {
        auto bCandidate = bSurfaceRep.intersect(gctx, position, sDirection,
                                                options.boundaryCheck);
        // Intersect and continue
        auto bIntersection = checkIntersection(bCandidate, bsIter.get());
        if (bIntersection) {
          ACTS_VERBOSE(" - Proceed with surface");
          bIntersections.push_back(bIntersection);
        } else {
          ACTS_VERBOSE(" - Surface intersecion invalid");
        }
      } else {
        ACTS_VERBOSE(" - Surface is excluded surface");
      }
    }
  };

  // Process the boundaries of the current volume
  auto& bSurfaces = boundarySurfaces();
  ACTS_VERBOSE("Volume reports " << bSurfaces.size() << " boundary surfaces");
  processBoundaries(bSurfaces);

  // Process potential boundaries of contained volumes
  auto confinedDenseVolumes = denseVolumes();
  ACTS_VERBOSE("Volume reports " << confinedDenseVolumes.size()
                                 << " confined dense volumes");
  for (const auto& dv : confinedDenseVolumes) {
    auto& bSurfacesConfined = dv->boundarySurfaces();
    ACTS_VERBOSE(" -> " << bSurfacesConfined.size() << " boundary surfaces");
    processBoundaries(bSurfacesConfined);
  }

  auto comparator = [](double a, double b) {
    // sign function would be nice but ...
    if ((a > 0 && b > 0) || (a < 0 && b < 0)) {
      return a < b;
    }
    if (a > 0) {  // b < 0
      return true;
    }
    return false;
  };

  // Sort them accordingly to the navigation direction
  if (options.navDir == NavigationDirection::Forward) {
    std::sort(bIntersections.begin(), bIntersections.end(),
              [&](const auto& a, const auto& b) {
                return comparator(a.intersection.pathLength,
                                  b.intersection.pathLength);
              });
  } else {
    std::sort(bIntersections.begin(), bIntersections.end(),
              [&](const auto& a, const auto& b) {
                return comparator(-a.intersection.pathLength,
                                  -b.intersection.pathLength);
              });
  }
  return bIntersections;
}

boost::container::small_vector<Acts::LayerIntersection, 10>
Acts::TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Layer>& options) const {
  // the layer intersections which are valid
  boost::container::small_vector<Acts::LayerIntersection, 10> lIntersections;

  // the confinedLayers
  if (m_confinedLayers != nullptr) {
    // start layer given or not - test layer
    const Layer* tLayer = options.startObject != nullptr
                              ? options.startObject
                              : associatedLayer(gctx, position);
    while (tLayer != nullptr) {
      // check if the layer needs resolving
      // - resolveSensitive -> always take layer if it has a surface array
      // - resolveMaterial -> always take layer if it has material
      // - resolvePassive -> always take, unless it's a navigation layer
      // skip the start object
      if (tLayer != options.startObject && tLayer->resolve(options)) {
        // if it's a resolveable start layer, you are by definition on it
        // layer on approach intersection
        auto atIntersection =
            tLayer->surfaceOnApproach(gctx, position, direction, options);
        auto path = atIntersection.intersection.pathLength;
        bool withinLimit = std::abs(path) <= std::abs(options.pathLimit);
        // Intersection is ok - take it (move to surface on appraoch)
        if (atIntersection &&
            (atIntersection.object != options.targetSurface) && withinLimit) {
          // create a layer intersection
          lIntersections.push_back(LayerIntersection(
              atIntersection.intersection, tLayer, atIntersection.object));
        }
      }
      // move to next one or break because you reached the end layer
      tLayer =
          (tLayer == options.endObject)
              ? nullptr
              : tLayer->nextLayer(gctx, position, options.navDir * direction);
    }
    // sort them accordingly to the navigation direction
    if (options.navDir == NavigationDirection::Forward) {
      std::sort(lIntersections.begin(), lIntersections.end());
    } else {
      std::sort(lIntersections.begin(), lIntersections.end(), std::greater<>());
    }
  }
  // and return
  return lIntersections;
}

namespace {
template <typename T>
std::vector<const Acts::Volume*> intersectSearchHierarchy(
    const T obj, const Acts::Volume::BoundingBox* lnode) {
  std::vector<const Acts::Volume*> hits;
  hits.reserve(20);  // arbitrary
  do {
    if (lnode->intersect(obj)) {
      if (lnode->hasEntity()) {
        // found primitive
        // check obb to limit false positivies
        const Acts::Volume* vol = lnode->entity();
        const auto& obb = vol->orientedBoundingBox();
        if (obb.intersect(obj.transformed(vol->itransform()))) {
          hits.push_back(vol);
        }
        // we skip in any case, whether we actually hit the OBB or not
        lnode = lnode->getSkip();
      } else {
        // go over children
        lnode = lnode->getLeftChild();
      }
    } else {
      lnode = lnode->getSkip();
    }
  } while (lnode != nullptr);

  return hits;
}
}  // namespace

std::vector<Acts::SurfaceIntersection>
Acts::TrackingVolume::compatibleSurfacesFromHierarchy(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, double angle,
    const NavigationOptions<Surface>& options) const {
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(20);  // arbitrary

  // The limits for this navigation step
  double pLimit = options.pathLimit;
  double oLimit = options.overstepLimit;

  if (m_bvhTop == nullptr) {
    return sIntersections;
  }

  // The signed direction
  Vector3 sdir = options.navDir * direction;

  std::vector<const Volume*> hits;
  if (angle == 0) {
    // use ray
    Ray3D obj(position, sdir);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  } else {
    Acts::Frustum<ActsScalar, 3, 4> obj(position, sdir, angle);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  }

  // have cells, decompose to surfaces
  for (const Volume* vol : hits) {
    const AbstractVolume* avol = dynamic_cast<const AbstractVolume*>(vol);
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        boundarySurfaces = avol->boundarySurfaces();
    for (const auto& bs : boundarySurfaces) {
      const Surface& srf = bs->surfaceRepresentation();
      auto sfi = srf.intersect(gctx, position, sdir, false);
      if (sfi and sfi.intersection.pathLength > oLimit and
          sfi.intersection.pathLength <= pLimit) {
        sIntersections.push_back(std::move(sfi));
      }
    }
  }

  // Sort according to the path length
  if (options.navDir == NavigationDirection::Forward) {
    std::sort(sIntersections.begin(), sIntersections.end());
  } else {
    std::sort(sIntersections.begin(), sIntersections.end(), std::greater<>());
  }

  return sIntersections;
}
