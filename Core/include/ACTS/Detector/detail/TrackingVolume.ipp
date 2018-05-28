// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.ipp, ACTS project
///////////////////////////////////////////////////////////////////

template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<LayerIntersection>
TrackingVolume::compatibleLayers(const parameters_t& parameters,
                                 const options_t&    options,
                                 const corrector_t&  corrfnc) const
{

  // get position and momentum from the parameters
  const Vector3D& pos = parameters.position();
  Vector3D        mom = options.navDir * parameters.momentum();

  // the layer intersections which are valid
  std::vector<LayerIntersection> lIntersections;
  // the confinedLayers
  if (m_confinedLayers) {
    // start layer given or not - test layer
    const Layer* tLayer
        = options.startObject ? options.startObject : associatedLayer(pos);
    if (tLayer) {
      do {
        // check if the layer needs resolving
        // - resolveSensitive -> always take layer if it has a surface array
        // - resolveMaterial -> always take layer if it has material
        // - resolvePassive -> always take, unless it's a navigation layer
        if (tLayer->resolve(options)) {
          // if it's a resolveable start layer, you ar by definition on it
          if (tLayer == options.startObject) {
            // create an intersection with path length 0.
            Intersection   cIntersection(pos, 0., true);
            const Surface* tSurface = &(tLayer->surfaceRepresentation());
            lIntersections.push_back(
                LayerIntersection(cIntersection, tLayer, tSurface));
          } else {
            // layer on approach intersection
            auto atIntersection
                = tLayer->surfaceOnApproach(parameters, options, corrfnc);
            auto path = atIntersection.intersection.pathLength;
            // Intersection is ok - take it (move to surface on appraoch)
            if (atIntersection
                && path * path <= options.pathLimit * options.pathLimit) {
              // create a layer intersection
              lIntersections.push_back(LayerIntersection(
                  atIntersection.intersection, tLayer, atIntersection.object));
            }
          }
        }
        // move to next one or break because you reached the end layer
        tLayer = (tLayer == options.endObject) ? nullptr
                                               : tLayer->nextLayer(pos, mom);
      } while (tLayer);
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
// straight line intersection @todo change hard-coded default
template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<BoundaryIntersection>
TrackingVolume::compatibleBoundaries(const parameters_t& parameters,
                                     const options_t&    options,
                                     const corrector_t&  corrfnc) const
{
  // Loop over boundarySurfaces and calculate the intersection
  std::vector<BoundaryIntersection> bIntersections;
  auto excludeObject = options.startObject;
  auto&                             bSurfaces = boundarySurfaces();
  for (auto& bsIter : bSurfaces) {
    // get the boundary surface pointer
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    const auto& bSurfaceRep = bSurface->surfaceRepresentation();
    // exclude the on boundary object
    if (excludeObject && excludeObject == &bSurfaceRep) continue;
    // intersect the surface 
    SurfaceIntersection bsIntersection
        = bSurfaceRep.intersectionEstimate(parameters, options, corrfnc);
    // check if the intersection is valid, but exlude the on-surface case
    // when requested -- move to intersectionestimate
    if (bsIntersection)
      bIntersections.push_back(BoundaryIntersection(
          bsIntersection.intersection, bSurface, &bSurfaceRep, options.navDir));
  }
  // and now sort to get the closest - need custom sort here to respect sign
  // sort them accordingly to the path length
  if (options.navDir == forward) {
    std::sort(bIntersections.begin(), bIntersections.end());
  } else {
    std::sort(bIntersections.begin(), bIntersections.end(), std::greater<>());
  }
  // and return
  return bIntersections;
}

template <typename parameters_t>
std::vector<LayerIntersection>
TrackingVolume::layerCandidatesOrdered(const Layer*         sLayer,
                                       const Layer*         eLayer,
                                       const parameters_t&  pars,
                                       NavigationDirection  pDir,
                                       const BoundaryCheck& bcheck,
                                       bool                 resolveSensitive,
                                       bool                 resolveMaterial,
                                       bool resolvePassive) const
{
  // get position and momentum from the parameters
  const Vector3D& gp = pars.position();
  const Vector3D& gm = pars.momentum();

  // the layer intersections
  std::vector<LayerIntersection> lIntersections;
  // assign the direction
  const Vector3D& dir
      = (pDir == forward) ? gm.unit() : Vector3D(-1 * gm.unit());
  // the confinedLayers
  if (m_confinedLayers) {
    // cache the longest path length to avoid punch-through to the other side
    Intersection   sLayerIntersection(Vector3D(0., 0., 0), 0., true, 0.);
    const Surface* sLayerSurface   = 0;
    double         validPathLength = 0.;
    // - get compatible layers back from the LayerArray simply because of the
    // binning
    // start layer given or not - test layer
    const Layer* tLayer = sLayer ? sLayer : associatedLayer(gp);
    if (tLayer) {
      do {
        // Fist:
        // - collect the startObject
        // Then:
        // - resolveSensitive -> always take layer if it has a surface array
        // - resolveMaterial -> always take layer if it has material
        // - resolvePassive -> always take, unless it's a navigation layer
        // Last:
        // - also collect the finalLayer
        if (tLayer->resolve(resolveSensitive, resolveMaterial, resolvePassive)
            || tLayer == sLayer
            || tLayer == eLayer) {
          // layer on approach intersection
          auto atIntersection = tLayer->surfaceOnApproach(gp,
                                                          gm,
                                                          pDir,
                                                          bcheck,
                                                          resolveSensitive,
                                                          resolveMaterial,
                                                          resolvePassive);

          // (a) if the current layer is NOT the start layer
          // - intersection is ok
          if (tLayer != sLayer && atIntersection.intersection.valid) {
            // create a layer intersection
            lIntersections.push_back(
                LayerIntersection(atIntersection.intersection,
                                  tLayer,
                                  atIntersection.object,
                                  pDir));
            validPathLength = atIntersection.intersection.pathLength;
          } else if (tLayer == sLayer) {
            // (b) the current layer is the start layer - we need to cache it
            // and check with the path length
            //     this avoids potential punch-through to other side of
            sLayerIntersection = atIntersection.intersection;
            sLayerSurface      = atIntersection.object;
          } else if (tLayer == eLayer) {
            // (c) it is the end layer after all
            // - provide it and break the loop
            lIntersections.push_back(
                LayerIntersection(atIntersection.intersection,
                                  tLayer,
                                  atIntersection.object,
                                  pDir));
            break;
          }
        }
        // move to next one or break because you reached the end layer
        tLayer = (tLayer == eLayer) ? nullptr : tLayer->nextLayer(gp, dir);
      } while (tLayer);
    }

    // final check for compatibility of the start layer in order to avoid
    // punch-through
    if (sLayer && sLayerIntersection.valid
        && sLayerIntersection.pathLength < validPathLength
        && sLayer->resolve(resolveSensitive, resolveMaterial, resolvePassive))
      lIntersections.push_back(
          LayerIntersection(sLayerIntersection, sLayer, sLayerSurface, pDir));

    // sort them accordingly to the path length
    std::sort(lIntersections.begin(), lIntersections.end());
  }
  // and return
  return lIntersections;
}

// Returns the boundary surfaces ordered in probability to hit them based on
// straight line intersection @todo change hard-coded default
template <typename parameters_t>
std::vector<BoundaryIntersection>
TrackingVolume::boundarySurfacesOrdered(const parameters_t& pars,
                                        NavigationDirection pDir,
                                        bool                skipCurrent) const
{
  // assign the direction
  const Vector3D mom
      = (pDir == forward ? pars.momentum() : Vector3D(-1 * pars.momentum()));
  // loop over boundarySurfaces and calculate the intersection
  std::vector<BoundaryIntersection> bIntersections;
  auto&                             bSurfaces = boundarySurfaces();
  for (auto& bsIter : bSurfaces) {
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    Intersection                            bsIntersection
        = bSurface->surfaceRepresentation().intersectionEstimate(
            pars.position(), mom, forward, false);
    // check if the intersection is valid, but exlude the on-surface case
    // when requested
    if (bsIntersection.valid
        && (!skipCurrent
            || std::abs(bsIntersection.pathLength) < s_onSurfaceTolerance))
      bIntersections.push_back(
          BoundaryIntersection(bsIntersection,
                               bSurface,
                               &(bSurface->surfaceRepresentation()),
                               pDir));
  }
  // and now sort to get the closest
  std::sort(bIntersections.begin(), bIntersections.end());
  // and return
  return bIntersections;
}
