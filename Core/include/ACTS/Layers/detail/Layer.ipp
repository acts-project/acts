// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <class T>
bool
Layer::onLayer(const T& pars, const BoundaryCheck& bcheck) const
{
  return isOnLayer(pars.position(), bcheck);
}

template <class T>
bool
Layer::getCompatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                             const T&                          pars,
                             PropDirection                     pDir,
                             const BoundaryCheck&              bcheck,
                             bool                              collectSensitive,
                             bool                              collectPassive,
                             int                               searchType,
                             const Surface*                    startSurface,
                             const Surface*                    endSurface,
                             const ICompatibilityEstimator*    ice) const
{
  // @TODO check if the approach surface should be provided, would be a waste
  // since it's already propagated for

  // intersectionTest with searchType
  bool intersectionTest = !(searchType % 2);

  // fast exit - nothing to do
  if (!m_surfaceArray || !m_approachDescriptor)
    return false;

  // position and momentum/dir
  const Vector3D& pos = pars.position();
  const Vector3D  dir = (pDir == oppositeMomentum)
      ? Vector3D(-1. * pars.momentum().unit())
      : pars.momentum().unit();

  // check if you have to stop at the endSurface
  double maxPathLength = 10e10;
  if (endSurface) {
    // intersect the end surface
    Intersection endInter
        = endSurface->intersectionEstimate(pos, dir, pDir, bcheck);
    // non-valid intersection with the end surface provided at this layer
    // indicates wrong direction or faulty setup
    // -> do not return compatible surfaces since they may lead you on a wrong
    // navigation path
    if (endInter.valid && endInter.pathLength > 0.)
      maxPathLength = endInter.pathLength;
    else
      return cSurfaces.size();
    // search Type will be increased automatically
    // when endSurface is given and no test done
    if (searchType % 2) ++searchType;
  }
  // we have a contained surface array
  // - resolve the different search modes
  // compatible test surfaces
  std::vector<const Surface*> ctestSurfaces;
  // this is providing all surfaces to the extrapolation eninge
  if (searchType <= 0) {
    // take all the test surfaces & their bin mates
    auto allTestSurfaces = m_surfaceArray->arrayObjects();
    // reserve twice the amount
    ctestSurfaces.reserve(allTestSurfaces.size());
    for (auto& atSurface : allTestSurfaces) {
      // get the bin mates if they exist
      if (atSurface && atSurface->associatedDetectorElement()) {
        // get the bin mates
        auto bmElements = atSurface->associatedDetectorElement()->binmembers();
        for (auto& bmElement : bmElements)
          ctestSurfaces.push_back(&(bmElement->surface()));
      }
      ctestSurfaces.push_back(atSurface);
    }
  } else if (m_surfaceArray) {
    // get the nominal test object
    auto tSurface = m_surfaceArray->object(pos);
    if (tSurface && tSurface->associatedDetectorElement()) {
      // get the detector elements
      auto dElement = tSurface->associatedDetectorElement();
      std::vector<const DetectorElementBase*> dElements = {dElement};
      // get the neighbours
      dElements.insert(dElements.begin(),
                       dElement->neighbours().begin(),
                       dElement->neighbours().end());
      // loop over all detector elements and add their surfaces and binmember
      // surfaces
      for (auto& ade : dElements) {
        // insert the surface
        ctestSurfaces.push_back(&(ade->surface()));
        // insert the bin members and the neighbors
        for (auto& abm : ade->binmembers()) {
          ctestSurfaces.push_back(&(abm->surface()));
        }
      }
    }
  }

  // loop over all the possible
  // sensitive surfaces and test them
  for (auto& ctSurface : ctestSurfaces)
    testCompatibleSurface(cSurfaces,
                          *ctSurface,
                          pos,
                          dir,
                          pDir,
                          bcheck,
                          maxPathLength,
                          collectSensitive,
                          collectPassive,
                          intersectionTest,
                          startSurface,
                          endSurface,
                          ice);

  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  testCompatibleSurface(cSurfaces,
                        *layerSurface,
                        pos,
                        dir,
                        pDir,
                        bcheck,
                        maxPathLength,
                        collectSensitive,
                        collectPassive,
                        intersectionTest,
                        startSurface,
                        endSurface,
                        ice);

  // the approach surfaces are always testSurfaces
  // usually, the surface on approach is excluded
  if (m_approachDescriptor) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces
        = m_approachDescriptor->containedSurfaces();
    for (auto& aSurface : approachSurfaces)
      testCompatibleSurface(cSurfaces,
                            *aSurface,
                            pos,
                            dir,
                            pDir,
                            bcheck,
                            maxPathLength,
                            collectSensitive,
                            collectPassive,
                            intersectionTest,
                            startSurface,
                            endSurface,
                            ice);
  }

  // only sort it if the intersection was done
  if (intersectionTest) std::sort(cSurfaces.begin(), cSurfaces.end());

  // return
  return intersectionTest;
}

inline void
Layer::testCompatibleSurface(std::vector<SurfaceIntersection>& cSurfaces,
                             const Surface&                    surface,
                             const Vector3D&                   pos,
                             const Vector3D&                   dir,
                             PropDirection                     pDir,
                             const BoundaryCheck&              bcheck,
                             double                            maxPathLength,
                             bool                              collectSensitive,
                             bool                              collectPassive,
                             bool                              intersectionTest,
                             const Surface*                    startSurface,
                             const Surface*                    endSurface,
                             const ICompatibilityEstimator*) const
{
  // fast exists
  // (1) skip the start and end surface
  if (&surface == endSurface || &surface == startSurface) return;
  // (2) no material, no passive collection, not active
  if (!surface.associatedDetectorElement() && !collectPassive
      && surface.associatedMaterial() == nullptr)
    return;
  // (3) is active, not configured for collect active
  if (surface.associatedDetectorElement() && !collectSensitive
      && surface.associatedMaterial() == nullptr)
    return;
  // then take it if you don't have to do an intersection
  if (!intersectionTest) {
    // return the surface
    cSurfaces.push_back(
        SurfaceIntersection(Intersection(pos, 0., true), &surface, pDir));
    // job done
    return;
  }
  // check if you need to force the momentum direction
  bool fDirection = (pDir == anyDirection ? false : true);
  // the intersection
  Intersection sfIntersection
      = surface.intersectionEstimate(pos, dir, fDirection, bcheck);
  // check if the intersection is valid and the maxPathLength has not been
  // exceeded
  if (sfIntersection.valid && sfIntersection.pathLength < maxPathLength) {
    // resulting propDirection
    PropDirection rDir
        = (sfIntersection.pathLength > 0 ? alongMomentum : oppositeMomentum);
    // and the surfaces & direction to push back - take all
    if (collectPassive
        || (collectSensitive && surface.associatedDetectorElement())
        || surface.associatedMaterial())
      cSurfaces.push_back(SurfaceIntersection(sfIntersection, &surface, rDir));
  }
}
