// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Acts::Layer* TrackingVolume::associatedLayer(
    const GeometryContext& /*gctx*/, const Vector3D& position) const {
  // confined static layers - highest hierarchy
  if (m_confinedLayers != nullptr) {
    return (m_confinedLayers->object(position).get());
  }

  // return the null pointer
  return nullptr;
}

template <typename options_t>
std::vector<LayerIntersection> TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const options_t& options) const {
  // the layer intersections which are valid
  std::vector<LayerIntersection> lIntersections;

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
        bool withinLimit =
            (path * path <= options.pathLimit * options.pathLimit);
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
    if (options.navDir == forward) {
      std::sort(lIntersections.begin(), lIntersections.end());
    } else {
      std::sort(lIntersections.begin(), lIntersections.end(), std::greater<>());
    }
  }
  // and return
  return lIntersections;
}

// Returns the boundary surfaces ordered in probability to hit them based on
template <typename options_t>
std::vector<BoundaryIntersection> TrackingVolume::compatibleBoundaries(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const options_t& options) const {
  // Loop over boundarySurfaces and calculate the intersection
  auto excludeObject = options.startObject;
  std::vector<BoundaryIntersection> bIntersections;

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

    double cLimit = sIntersection.intersection.pathLength;
    // Check if the surface is within limit
    bool withinLimit =
        (cLimit > oLimit and
         cLimit * cLimit <= pLimit * pLimit + s_onSurfaceTolerance);
    if (withinLimit) {
      sIntersection.intersection.pathLength *=
          std::copysign(1., options.navDir);
      return BoundaryIntersection(sIntersection.intersection, bSurface,
                                  sIntersection.object);
    }
    // Check the alternative
    if (sIntersection.alternative) {
      // Test the alternative
      cLimit = sIntersection.alternative.pathLength;
      withinLimit = (cLimit > oLimit and
                     cLimit * cLimit <= pLimit * pLimit + s_onSurfaceTolerance);
      if (sIntersection.alternative and withinLimit) {
        sIntersection.alternative.pathLength *=
            std::copysign(1., options.navDir);
        return BoundaryIntersection(sIntersection.alternative, bSurface,
                                    sIntersection.object);
      }
    }
    // Return an invalid one
    return BoundaryIntersection();
  };

  /// Helper function to process boundary surfaces
  auto processBoundaries =
      [&](const TrackingVolumeBoundaries& bSurfaces) -> void {
    // Loop over the boundary surfaces
    for (auto& bsIter : bSurfaces) {
      // Get the boundary surface pointer
      const auto& bSurfaceRep = bsIter->surfaceRepresentation();
      // Exclude the boundary where you are on
      if (excludeObject != &bSurfaceRep) {
        auto bCandidate = bSurfaceRep.intersect(gctx, position, sDirection,
                                                options.boundaryCheck);
        // Intersect and continue
        auto bIntersection = checkIntersection(bCandidate, bsIter.get());
        if (bIntersection) {
          bIntersections.push_back(bIntersection);
        }
      }
    }
  };

  // Process the boundaries of the current volume
  auto& bSurfaces = boundarySurfaces();
  processBoundaries(bSurfaces);

  // Process potential boundaries of contained volumes
  auto confinedDenseVolumes = denseVolumes();
  for (const auto& dv : confinedDenseVolumes) {
    auto& bSurfacesConfined = dv->boundarySurfaces();
    processBoundaries(bSurfacesConfined);
  }

  // Sort them accordingly to the navigation direction
  if (options.navDir == forward) {
    std::sort(bIntersections.begin(), bIntersections.end());
  } else {
    std::sort(bIntersections.begin(), bIntersections.end(), std::greater<>());
  }
  return bIntersections;
}

template <typename options_t>
std::vector<SurfaceIntersection>
TrackingVolume::compatibleSurfacesFromHierarchy(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, double angle, const options_t& options) const {
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(20);  // arbitrary

  // The limits for this navigation step
  double pLimit = options.pathLimit;
  double oLimit = options.overstepLimit;

  if (m_bvhTop == nullptr || !options.navDir) {
    return sIntersections;
  }

  // The signed direction
  Vector3D sdir = options.navDir * direction;

  std::vector<const Volume*> hits;
  if (angle == 0) {
    // use ray
    Ray3D obj(position, sdir);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  } else {
    Acts::Frustum<double, 3, 4> obj(position, sdir, angle);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  }

  // have cells, decompose to surfaces
  for (const Volume* vol : hits) {
    const AbstractVolume* avol = dynamic_cast<const AbstractVolume*>(vol);
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        boundarySurfaces = avol->boundarySurfaces();
    for (const auto& bs : boundarySurfaces) {
      const Surface& srf = bs->surfaceRepresentation();
      SurfaceIntersection sfi(
          srf.intersectionEstimate(gctx, position, sdir, false), &srf);

      if (sfi and sfi.intersection.pathLength > oLimit and
          sfi.intersection.pathLength <= pLimit) {
        sIntersections.push_back(std::move(sfi));
      }
    }
  }

  // Sort according to the path length
  if (options.navDir == forward) {
    std::sort(sIntersections.begin(), sIntersections.end());
  } else {
    std::sort(sIntersections.begin(), sIntersections.end(), std::greater<>());
  }

  return sIntersections;
}

template <typename T>
std::vector<const Volume*> TrackingVolume::intersectSearchHierarchy(
    const T obj, const Volume::BoundingBox* lnode) {
  std::vector<const Volume*> hits;
  hits.reserve(20);  // arbitrary
  do {
    if (lnode->intersect(obj)) {
      if (lnode->hasEntity()) {
        // found primitive
        // check obb to limit false positivies
        const Volume* vol = lnode->entity();
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
