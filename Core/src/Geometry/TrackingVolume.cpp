// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingVolume.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GlueVolumesDescriptor.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Ray.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <ostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace Acts {
class ISurfaceMaterial;

TrackingVolume::TrackingVolume(
    const Transform3& transform, std::shared_ptr<VolumeBounds> volbounds,
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
TrackingVolume::TrackingVolume(
    const Transform3& transform, std::shared_ptr<VolumeBounds> volumeBounds,
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

TrackingVolume::~TrackingVolume() {
  delete m_glueVolumeDescriptor;
}

const TrackingVolume* TrackingVolume::lowestTrackingVolume(
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

const TrackingVolumeBoundaries& TrackingVolume::boundarySurfaces() const {
  return (m_boundarySurfaces);
}

void TrackingVolume::connectDenseBoundarySurfaces(
    MutableTrackingVolumeVector& confinedDenseVolumes) {
  if (!confinedDenseVolumes.empty()) {
    Direction dir = Direction::Positive;
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
          dir = Direction::Positive;
          mutableBs->attachVolume(this, dir);
        } else {
          if (mutableBs->m_oppositeVolume == nullptr &&
              mutableBs->m_alongVolume != nullptr) {
            dir = Direction::Negative;
            mutableBs->attachVolume(this, dir);
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

void TrackingVolume::createBoundarySurfaces() {
  using Boundary = BoundarySurfaceT<TrackingVolume>;

  // Transform Surfaces To BoundarySurfaces
  auto orientedSurfaces = Volume::volumeBounds().orientedSurfaces(m_transform);

  m_boundarySurfaces.reserve(orientedSurfaces.size());
  for (auto& osf : orientedSurfaces) {
    TrackingVolume* opposite = nullptr;
    TrackingVolume* along = nullptr;
    if (osf.direction == Direction::OppositeNormal) {
      opposite = this;
    } else {
      along = this;
    }
    m_boundarySurfaces.push_back(std::make_shared<const Boundary>(
        std::move(osf.surface), opposite, along));
  }
}

void TrackingVolume::glueTrackingVolume(const GeometryContext& gctx,
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
  // Coerce the arbitrary position bPosition to be on the surface repr so we can
  // get a normal
  Vector3 nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  Direction dir = Direction::fromScalar(nvector.dot(distance));
  // The easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr) ||
      m_glueVolumeDescriptor->glueVolumes(bsfMine) == nullptr) {
    // the boundary orientation
    auto mutableBSurfaceMine =
        std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(bSurfaceMine);
    mutableBSurfaceMine->attachVolume(neighbor, dir);
    // Make sure you keep the boundary material if there
    const Surface& neighborSurface =
        neighbor->m_boundarySurfaces.at(bsfNeighbor)->surfaceRepresentation();
    auto neighborMaterial = neighborSurface.surfaceMaterialSharedPtr();
    const Surface& mySurface = bSurfaceMine->surfaceRepresentation();
    auto myMaterial = mySurface.surfaceMaterialSharedPtr();
    // Keep the neighbor material
    if (myMaterial == nullptr && neighborMaterial != nullptr) {
      Surface* myMutbableSurface = const_cast<Surface*>(&mySurface);
      myMutbableSurface->assignSurfaceMaterial(neighborMaterial);
    }
    // Now set it to the neighbor volume
    (neighbor->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
  }
}

void TrackingVolume::glueTrackingVolumes(
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
  // Coerce the arbitrary position bPosition to be on the surface repr so we can
  // get a normal
  Vector3 nvector =
      bSurfaceMine->surfaceRepresentation().normal(gctx, bPosition);
  // estimate the orientation
  Direction dir = Direction::fromScalar(nvector.dot(distance));
  // the easy case :
  // - no glue volume descriptors on either side
  if ((m_glueVolumeDescriptor == nullptr) ||
      !m_glueVolumeDescriptor->glueVolumes(bsfMine)) {
    // the boundary orientation
    auto mutableBSurfaceMine =
        std::const_pointer_cast<BoundarySurfaceT<TrackingVolume>>(bSurfaceMine);
    mutableBSurfaceMine->attachVolumeArray(neighbors, dir);
    // now set it to the neighbor volumes - the optised way
    for (auto& nVolume : neighbors->arrayObjects()) {
      auto mutableNVolume = std::const_pointer_cast<TrackingVolume>(nVolume);
      (mutableNVolume->m_boundarySurfaces).at(bsfNeighbor) = bSurfaceMine;
    }
  }
}

void TrackingVolume::assignBoundaryMaterial(
    std::shared_ptr<const ISurfaceMaterial> surfaceMaterial,
    BoundarySurfaceFace bsFace) {
  auto bSurface = m_boundarySurfaces.at(bsFace);
  RegularSurface* surface =
      const_cast<RegularSurface*>(&bSurface->surfaceRepresentation());
  surface->assignSurfaceMaterial(std::move(surfaceMaterial));
}

void TrackingVolume::updateBoundarySurface(
    BoundarySurfaceFace bsf,
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs,
    bool checkmaterial) {
  if (checkmaterial) {
    auto cMaterialPtr = m_boundarySurfaces.at(bsf)
                            ->surfaceRepresentation()
                            .surfaceMaterialSharedPtr();
    auto bsMaterial = bs->surfaceRepresentation().surfaceMaterial();
    if (cMaterialPtr != nullptr && bsMaterial == nullptr) {
      RegularSurface* surface =
          const_cast<RegularSurface*>(&bs->surfaceRepresentation());
      surface->assignSurfaceMaterial(cMaterialPtr);
    }
  }
  m_boundarySurfaces.at(bsf) = std::move(bs);
}

void TrackingVolume::registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd) {
  delete m_glueVolumeDescriptor;
  m_glueVolumeDescriptor = gvd;
}

GlueVolumesDescriptor& TrackingVolume::glueVolumesDescriptor() {
  if (m_glueVolumeDescriptor == nullptr) {
    m_glueVolumeDescriptor = new GlueVolumesDescriptor;
  }
  return (*m_glueVolumeDescriptor);
}

void TrackingVolume::synchronizeLayers(double envelope) const {
  // case a : Layers exist
  // msgstream << MSG::VERBOSE << "  -> synchronizing Layer dimensions of
  // TrackingVolume '" << volumeName() << "'." << endreq;

  if (m_confinedLayers) {
    // msgstream << MSG::VERBOSE << "  ---> working on " <<
    // m_confinedLayers->arrayObjects().size() << " (material+navigation)
    // layers." << endreq;
    for (auto& clayIter : m_confinedLayers->arrayObjects()) {
      if (clayIter) {
        // @todo implement synchronize layer
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

void TrackingVolume::interlinkLayers() {
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

void TrackingVolume::closeGeometry(
    const IMaterialDecorator* materialDecorator,
    std::unordered_map<GeometryIdentifier, const TrackingVolume*>& volumeMap,
    std::size_t& vol, const GeometryIdentifierHook& hook,
    const Logger& logger) {
  // we can construct the volume ID from this
  auto volumeID = GeometryIdentifier().setVolume(++vol);
  // assign the Volume ID to the volume itself
  auto thisVolume = const_cast<TrackingVolume*>(this);
  thisVolume->assignGeometryId(volumeID);
  ACTS_DEBUG("volumeID: " << volumeID << ", name: " << volumeName());
  // insert the volume into the map
  volumeMap[volumeID] = thisVolume;

  // assign the material if you have a decorator
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*thisVolume);
  }
  if (thisVolume->volumeMaterial() == nullptr &&
      thisVolume->motherVolume() != nullptr &&
      thisVolume->motherVolume()->volumeMaterial() != nullptr) {
    auto protoMaterial = dynamic_cast<const ProtoVolumeMaterial*>(
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
    // get the intersection solution
    auto& bSurface = bSurfIter->surfaceRepresentation();
    // create the boundary surface id
    auto boundaryID = GeometryIdentifier(volumeID).setBoundary(++iboundary);
    // now assign to the boundary surface
    auto& mutableBSurface = *(const_cast<RegularSurface*>(&bSurface));
    mutableBSurface.assignGeometryId(boundaryID);
    // Assign material if you have a decorator
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
        mutableLayerPtr->closeGeometry(materialDecorator, layerID, hook,
                                       logger);
      }
    }
  } else {
    // B) this is a container volume, go through sub volume
    // do the loop
    for (auto& volumesIter : m_confinedVolumes->arrayObjects()) {
      auto mutableVolumesIter =
          std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->setMotherVolume(this);
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol, hook,
                                        logger);
    }
  }

  if (!m_confinedDenseVolumes.empty()) {
    for (auto& volumesIter : m_confinedDenseVolumes) {
      auto mutableVolumesIter =
          std::const_pointer_cast<TrackingVolume>(volumesIter);
      mutableVolumesIter->setMotherVolume(this);
      mutableVolumesIter->closeGeometry(materialDecorator, volumeMap, vol, hook,
                                        logger);
    }
  }
}

// Returns the boundary surfaces ordered in probability to hit them based on
boost::container::small_vector<BoundaryIntersection, 4>
TrackingVolume::compatibleBoundaries(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction,
                                     const NavigationOptions<Surface>& options,
                                     const Logger& logger) const {
  ACTS_VERBOSE("Finding compatibleBoundaries");
  // Loop over boundarySurfaces and calculate the intersection
  auto excludeObject = options.startObject;
  boost::container::small_vector<BoundaryIntersection, 4> bIntersections;

  // The Limits: current, path & overstepping
  double nearLimit = 0;
  double farLimit = options.farLimit;

  // Helper function to test intersection
  auto checkIntersection =
      [&](SurfaceMultiIntersection& smIntersection,
          const BoundarySurface* bSurface) -> BoundaryIntersection {
    for (const auto& sIntersection : smIntersection.split()) {
      if (!sIntersection) {
        continue;
      }

      if (options.forceIntersectBoundaries) {
        const bool coCriterion =
            std::abs(sIntersection.pathLength()) < std::abs(nearLimit);
        ACTS_VERBOSE("Forcing intersection with surface "
                     << bSurface->surfaceRepresentation().geometryId());
        if (coCriterion) {
          ACTS_VERBOSE("Intersection forced successfully ");
          ACTS_VERBOSE("- intersection path length "
                       << std::abs(sIntersection.pathLength())
                       << " < overstep limit " << std::abs(nearLimit));
          return BoundaryIntersection(sIntersection, bSurface);
        }
        ACTS_VERBOSE("Can't force intersection: ");
        ACTS_VERBOSE("- intersection path length "
                     << std::abs(sIntersection.pathLength())
                     << " > overstep limit " << std::abs(nearLimit));
      }

      ACTS_VERBOSE("Check intersection with surface "
                   << bSurface->surfaceRepresentation().geometryId());
      if (detail::checkIntersection(sIntersection.intersection(), nearLimit,
                                    farLimit, logger)) {
        return BoundaryIntersection(sIntersection, bSurface);
      }
    }

    ACTS_VERBOSE("No intersection accepted");
    return BoundaryIntersection(SurfaceIntersection::invalid(), bSurface);
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
        auto bCandidate = bSurfaceRep.intersect(gctx, position, direction,
                                                options.boundaryCheck);
        // Intersect and continue
        auto bIntersection = checkIntersection(bCandidate, bsIter.get());
        if (bIntersection.first) {
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

  std::sort(bIntersections.begin(), bIntersections.end(),
            [&](const BoundaryIntersection& a, const BoundaryIntersection& b) {
              return comparator(a.first.pathLength(), b.first.pathLength());
            });
  return bIntersections;
}

boost::container::small_vector<LayerIntersection, 10>
TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Layer>& options) const {
  // the layer intersections which are valid
  boost::container::small_vector<LayerIntersection, 10> lIntersections;

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
        auto path = atIntersection.pathLength();
        bool withinLimit = std::abs(path) <= std::abs(options.farLimit);
        // Intersection is ok - take it (move to surface on approach)
        if (atIntersection && withinLimit) {
          // create a layer intersection
          lIntersections.push_back(LayerIntersection(atIntersection, tLayer));
        }
      }
      // move to next one or break because you reached the end layer
      tLayer = (tLayer == options.endObject)
                   ? nullptr
                   : tLayer->nextLayer(gctx, position, direction);
    }
    std::sort(lIntersections.begin(), lIntersections.end(),
              [](const LayerIntersection& a, const LayerIntersection& b) {
                return SurfaceIntersection::pathLengthOrder(a.first, b.first);
              });
  }
  // and return
  return lIntersections;
}

const std::string& TrackingVolume::volumeName() const {
  return m_name;
}

const IVolumeMaterial* TrackingVolume::volumeMaterial() const {
  return m_volumeMaterial.get();
}

const std::shared_ptr<const IVolumeMaterial>&
TrackingVolume::volumeMaterialSharedPtr() const {
  return m_volumeMaterial;
}

void TrackingVolume::assignVolumeMaterial(
    std::shared_ptr<const IVolumeMaterial> material) {
  m_volumeMaterial = std::move(material);
}

const LayerArray* TrackingVolume::confinedLayers() const {
  return m_confinedLayers.get();
}

const MutableTrackingVolumeVector TrackingVolume::denseVolumes() const {
  return m_confinedDenseVolumes;
}

std::shared_ptr<const TrackingVolumeArray> TrackingVolume::confinedVolumes()
    const {
  return m_confinedVolumes;
}

const TrackingVolume* TrackingVolume::motherVolume() const {
  return m_motherVolume;
}

TrackingVolume* TrackingVolume::motherVolume() {
  return m_motherVolume;
}

void TrackingVolume::setMotherVolume(TrackingVolume* mvol) {
  m_motherVolume = mvol;
}

const Acts::Layer* TrackingVolume::associatedLayer(
    const GeometryContext& /*gctx*/, const Vector3& position) const {
  // confined static layers - highest hierarchy
  if (m_confinedLayers != nullptr) {
    return (m_confinedLayers->object(position).get());
  }

  // return the null pointer
  return nullptr;
}
}  // namespace Acts
