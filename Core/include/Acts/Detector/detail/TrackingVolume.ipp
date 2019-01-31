// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.ipp, Acts project
///////////////////////////////////////////////////////////////////

template <typename options_t, typename corrector_t>
std::vector<LayerIntersection>
TrackingVolume::compatibleLayers(const Vector3D&    position,
                                 const Vector3D&    direction,
                                 const options_t&   options,
                                 const corrector_t& corrfnc) const
{

  // the layer intersections which are valid
  std::vector<LayerIntersection> lIntersections;

  // the confinedLayers
  if (m_confinedLayers) {
    // start layer given or not - test layer
    const Layer* tLayer
        = options.startObject ? options.startObject : associatedLayer(position);
    while (tLayer != nullptr) {
      // check if the layer needs resolving
      // - resolveSensitive -> always take layer if it has a surface array
      // - resolveMaterial -> always take layer if it has material
      // - resolvePassive -> always take, unless it's a navigation layer
      // skip the start object
      if (tLayer != options.startObject && tLayer->resolve(options)) {
        // if it's a resolveable start layer, you are by definition on it
        // layer on approach intersection
        auto atIntersection
            = tLayer->surfaceOnApproach(position, direction, options, corrfnc);
        auto path = atIntersection.intersection.pathLength;
        bool withinLimit
            = (path * path <= options.pathLimit * options.pathLimit);
        // Intersection is ok - take it (move to surface on appraoch)
        if (atIntersection && (atIntersection.object != options.targetSurface)
            && withinLimit) {
          // create a layer intersection
          lIntersections.push_back(LayerIntersection(
              atIntersection.intersection, tLayer, atIntersection.object));
        }
      }
      // move to next one or break because you reached the end layer
      tLayer = (tLayer == options.endObject)
          ? nullptr
          : tLayer->nextLayer(position, options.navDir * direction);
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
std::vector<LayerIntersection>
TrackingVolume::compatibleLayers(const parameters_t& parameters,
                                 const options_t&    options,
                                 const corrector_t&  corrfnc) const
{
  return compatibleLayers(
      parameters.position(), parameters.direction(), options, corrfnc);
}

// Returns the boundary surfaces ordered in probability to hit them based on
// straight line intersection @todo change hard-coded default
template <typename options_t, typename corrector_t, typename sorter_t>
std::vector<BoundaryIntersection>
TrackingVolume::compatibleBoundaries(const Vector3D&    position,
                                     const Vector3D&    direction,
                                     const options_t&   options,
                                     const corrector_t& corrfnc,
                                     const sorter_t&    sorter) const
{
  // Loop over boundarySurfaces and calculate the intersection
  auto  excludeObject = options.startObject;
  auto& bSurfaces     = boundarySurfaces();
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
  return sorter(nonExcludedBoundaries, position, direction, options, corrfnc);
}

// Returns the boundary surfaces ordered in probability to hit them based on
// straight line intersection @todo change hard-coded default
template <typename parameters_t,
          typename options_t,
          typename corrector_t,
          typename sorter_t>
std::vector<BoundaryIntersection>
TrackingVolume::compatibleBoundaries(const parameters_t& parameters,
                                     const options_t&    options,
                                     const corrector_t&  corrfnc,
                                     const sorter_t&     sorter) const
{
  return compatibleBoundaries(
      parameters.position(), parameters.direction(), options, corrfnc, sorter);
}