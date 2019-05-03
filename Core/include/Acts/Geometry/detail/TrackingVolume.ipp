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
    const GeometryContext& /*gctx*/, const Vector3D& gp) const {
  // confined static layers - highest hierarchy
  if (m_confinedLayers) {
    return (m_confinedLayers->object(gp).get());
  }

  // return the null pointer
  return nullptr;
}

template <typename options_t, typename corrector_t>
std::vector<LayerIntersection> TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const options_t& options,
    const corrector_t& corrfnc) const {
  // the layer intersections which are valid
  std::vector<LayerIntersection> lIntersections;

  // the confinedLayers
  if (m_confinedLayers) {
    // start layer given or not - test layer
    const Layer* tLayer = options.startObject ? options.startObject
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
        auto atIntersection = tLayer->surfaceOnApproach(
            gctx, position, direction, options, corrfnc);
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

template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<LayerIntersection> TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const parameters_t& parameters,
    const options_t& options, const corrector_t& corrfnc) const {
  return compatibleLayers(gctx, parameters.position(), parameters.direction(),
                          options, corrfnc);
}

// Returns the boundary surfaces ordered in probability to hit them based on
template <typename options_t, typename corrector_t, typename sorter_t>
std::vector<BoundaryIntersection> TrackingVolume::compatibleBoundaries(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const options_t& options,
    const corrector_t& corrfnc, const sorter_t& sorter) const {
  // Loop over boundarySurfaces and calculate the intersection
  auto excludeObject = options.startObject;
  auto& bSurfaces = boundarySurfaces();
  std::vector<const BoundarySurfaceT<TrackingVolume>*> nonExcludedBoundaries;

  for (auto& bsIter : bSurfaces) {
    // get the boundary surface pointer
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    const auto& bSurfaceRep = bSurface->surfaceRepresentation();
    // exclude the on boundary object
    if (excludeObject && excludeObject == &bSurfaceRep) {
      continue;
    }
    nonExcludedBoundaries.push_back(bSurface);
  }

  const std::vector<std::shared_ptr<TrackingVolume>> confinedDenseVolumes =
      denseVolumes();
  for (const auto& dv : confinedDenseVolumes) {
    auto& bSurfacesConfined = dv->boundarySurfaces();
    for (auto& bsIter : bSurfacesConfined) {
      // get the boundary surface pointer
      const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
      const auto& bSurfaceRep = bSurface->surfaceRepresentation();
      // exclude the on boundary object
      if (excludeObject && excludeObject == &bSurfaceRep) {
        continue;
      }
      nonExcludedBoundaries.push_back(bSurface);
    }
  }

  return sorter(gctx, nonExcludedBoundaries, position, direction, options,
                corrfnc);
}

// Returns the boundary surfaces ordered in probability to hit them based on
// straight line intersection @todo change hard-coded default
template <typename parameters_t, typename options_t, typename corrector_t,
          typename sorter_t>
std::vector<BoundaryIntersection> TrackingVolume::compatibleBoundaries(
    const GeometryContext& gctx, const parameters_t& parameters,
    const options_t& options, const corrector_t& corrfnc,
    const sorter_t& sorter) const {
  return compatibleBoundaries(gctx, parameters.position(),
                              parameters.direction(), options, corrfnc, sorter);
}

template <typename options_t, typename corrector_t>
std::vector<SurfaceIntersection>
TrackingVolume::compatibleSurfacesFromHierarchy(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, double angle, const options_t& options,
    const corrector_t& corrfnc) const {
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(20);  // arbitrary

  if (m_bvhTop == nullptr || !options.navDir) {
    return sIntersections;
  }

  Vector3D dir = direction;
  if (options.navDir == backward) {
    dir *= -1;
  }

  std::vector<const Volume*> hits;
  if (angle == 0) {
    // use ray
    Ray3D obj(position, dir);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  } else {
    Acts::Frustum<double, 3, 4> obj(position, dir, angle);
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
          srf.intersectionEstimate(gctx, position, direction, options.navDir,
                                   false, corrfnc),
          &srf);

      if (sfi) {
        sIntersections.push_back(std::move(sfi));
      }
    }
  }

  // sort according to the path length
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
