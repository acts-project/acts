// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>

namespace Acts {

inline const SurfaceArray*
Layer::surfaceArray() const
{
  return m_surfaceArray.get();
}

inline SurfaceArray*
Layer::surfaceArray()
{
  return const_cast<SurfaceArray*>(m_surfaceArray.get());
}

inline double
Layer::thickness() const
{
  return m_layerThickness;
}

inline LayerType
Layer::layerType() const
{
  return m_layerType;
}

inline const TrackingVolume*
Layer::trackingVolume() const
{
  return m_trackingVolume;
}

inline void
Layer::encloseTrackingVolume(const TrackingVolume& tvol)
{
  m_trackingVolume = &(tvol);
}

inline const DetachedTrackingVolume*
Layer::enclosingDetachedTrackingVolume() const
{
  return m_enclosingDetachedTrackingVolume;
}

inline void
Layer::encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol)
{
  m_enclosingDetachedTrackingVolume = &(tvol);
}

inline const AbstractVolume*
Layer::representingVolume() const
{
  return m_representingVolume;
}

inline const Layer*
Layer::nextLayer(const Vector3D& gp, const Vector3D& mom) const
{
  // no binutility -> no chance to find out the direction
  if (!m_nextLayerUtility) return nullptr;
  return (m_nextLayerUtility->nextDirection(gp, mom) < 0) ? m_nextLayers.first
                                                          : m_nextLayers.second;
}

inline bool
Layer::resolve(bool resolveSensitive,
               bool resolveMaterial,
               bool resolvePassive) const
{
  if (resolvePassive) return true;
  if (resolveSensitive && m_surfaceArray) return true;
  if (resolveMaterial && (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1
                          || surfaceRepresentation().associatedMaterial()))
    return true;
  return false;
}

inline void
Layer::registerRepresentingVolume(const AbstractVolume* theVol)
{
  delete m_representingVolume;
  m_representingVolume = theVol;
}

template <typename parameters_t>
bool
Layer::onLayer(const parameters_t& pars, const BoundaryCheck& bcheck) const
{
  return isOnLayer(pars.position(), bcheck);
}

template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<SurfaceIntersection>
Layer::compatibleSurfaces(const parameters_t& parameters,
                          const options_t&    options,
                          const corrector_t&  corrfnc) const
{
  // the list of valid intersection
  std::vector<SurfaceIntersection> sIntersections;
  // remember the surfaces for duplicate removal
  std::map<const Surface*, bool> accepted;

  // fast exit - there is nothing to
  if (!m_surfaceArray || !m_approachDescriptor || !options.navDir)
    return sIntersections;

  // reserve a few bins
  sIntersections.reserve(20);

  // (0) End surface check
  // @todo: - we might be able to skip this by use of options.pathLimit
  // check if you have to stop at the endSurface
  double maxPath = options.pathLimit;
  if (options.endObject) {
    // intersect the end surface
    // - it is the final one don't use the bounday check at all
    SurfaceIntersection endInter
        = options.endObject->template intersectionEstimate(
            parameters, options, corrfnc);
    // non-valid intersection with the end surface provided at this layer
    // indicates wrong direction or faulty setup
    // -> do not return compatible surfaces since they may lead you on a wrong
    // navigation path
    if (endInter)
      maxPath = endInter.intersection.pathLength;
    else
      return sIntersections;

  } else {
    // compatibleSurfaces() should only be called when on the layer,
    // i.e. the maximum path limit is given by the layer thickness times
    // path correction, we take a safety factor of 1.5
    // -> this avoids punch through for cylinders
    double pCorrection = surfaceRepresentation().pathCorrection(
        parameters.position(), parameters.momentum());
    maxPath = 1.5 * thickness() * pCorrection * options.navDir;
  }

  // lemma 0 : accept the surface
  auto acceptSurface = [&options, &accepted](const Surface& sf,
                                             bool sensitive = false) -> bool {
    // check for duplicates
    if (accepted.find(&sf) != accepted.end()) return false;
    // surface is sensitive and you're asked to resolve
    if (sensitive && options.resolveSensitive) return true;
    // next option: it's a material surface and you want to have it
    if (options.resolveMaterial && sf.associatedMaterial()) return true;
    // last option: resovle all
    return options.resolvePassive;
  };

  // lemma 1 : check and fill the surface
  // [&sIntersections, &options, &parameters,&corrfnc
  auto processSurface = [&](const Surface& sf, bool sensitive = false) {
    // veto if it's start or end surface
    if (options.startObject == &sf || options.endObject == &sf) return;
    // veto if it doesn't fit the prescription
    if (!acceptSurface(sf, sensitive)) return;
    // the surface intersection
    SurfaceIntersection sfi
        = sf.intersectionEstimate(parameters, options, corrfnc);
    // check if intersection is valid and pathLimit has not been exceeded
    double sifPath = sfi.intersection.pathLength;
    // check the maximum path length
    if (sfi && sifPath * sifPath <= maxPath * maxPath) {
      sIntersections.push_back(sfi);
      accepted[&sf] = true;
    }
    return;
  };

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  if (m_approachDescriptor
      && (options.resolveMaterial || options.resolvePassive)) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces
        = m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& aSurface : approachSurfaces) {
      processSurface(*aSurface);
    }
  }

  // (B) sensitive surface section
  //
  // check the sensitive surfaces if you have some
  if (m_surfaceArray && (options.resolveMaterial || options.resolvePassive
                         || options.resolveSensitive)) {
    // get the canditates
    const std::vector<const Surface*>& sensitiveSurfaces
        = m_surfaceArray->neighbors(parameters.position());
    // loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& sSurface : sensitiveSurfaces) {
      processSurface(*sSurface);
    }
  }

  // (C) representing surface section
  //
  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  processSurface(*layerSurface);

  // sort according to the path length
  if (options.navDir == forward)
    std::sort(sIntersections.begin(), sIntersections.end());
  else
    std::sort(sIntersections.begin(), sIntersections.end(), std::greater<>());

  return sIntersections;
}

template <typename parameters_t, typename options_t, typename corrector_t>
const SurfaceIntersection
Layer::surfaceOnApproach(const parameters_t& parameters,
                         const options_t&    options,
                         const corrector_t&  corrfnc) const
{
  // resolve directive based by options
  // - options.resolvePassive is on -> always
  // - options.resolveSensitive is on -> always
  // - options.resolveMaterial is on
  //   && either sensitive or approach surfaces have material
  bool resolvePS = options.resolveSensitive || options.resolvePassive;
  bool resolveMS = options.resolveMaterial
      && (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1
          || surfaceRepresentation().associatedMaterial());

  // now of course this only counts when you have an approach descriptor
  if (m_approachDescriptor && (resolvePS || resolveMS)) {
    // test if you are on an approach surface already, if so - provide it
    for (auto& asf : m_approachDescriptor->containedSurfaces()) {
      // in a connected geometry this is only a pointer comparison
      if (options.startObject
          && asf == &(options.startObject->surfaceRepresentation())) {
        Intersection nIntersection(parameters.position(), 0., true);
        return SurfaceIntersection(nIntersection, asf, options.navDir);
      }
    }
    // that's the collect trigger for always collecting
    // let's find the most suitable approach surface
    SurfaceIntersection aSurface
        = m_approachDescriptor->approachSurface(parameters, options, corrfnc);
    if (aSurface.intersection.valid) return (aSurface);
  }

  const Surface& rSurface = surfaceRepresentation();

  // if we have no approach descriptor - we have no sensitive surfaces
  if (rSurface.onSurface(parameters, options.boundaryCheck)) {
    Intersection nIntersection(parameters.position(), 0., true);
    return SurfaceIntersection(nIntersection, &rSurface, options.navDir);
  }

  // create the intersection with the surface representation
  return rSurface.intersectionEstimate(parameters, options, corrfnc);
}

// ----------- legacy method block: start ----------------------

template <typename parameters_t>
std::vector<SurfaceIntersection>
Layer::getCompatibleSurfaces(const parameters_t&  pars,
                             NavigationDirection  pDir,
                             const BoundaryCheck& bcheck,
                             bool                 resolveSensitive,
                             bool                 resolveMaterial,
                             bool                 resolvePassive,
                             int                  searchType,
                             const Surface*       startSurface,
                             const Surface*       endSurface) const
{

  // prepare the surface intersections for return
  std::vector<SurfaceIntersection> cSurfaces;
  cSurfaces.reserve(20);

  // fast exit - there is nothing to do nothing to do
  if (!m_surfaceArray || !m_approachDescriptor) return cSurfaces;

  // the test boudnary check is defined by the search type
  BoundaryCheck tCheck = searchType % 2 ? BoundaryCheck(false) : bcheck;

  // position and momentum/dir
  const Vector3D& pos = pars.position();
  const Vector3D  dir = (pDir == backward)
      ? Vector3D(-1. * pars.momentum().unit())
      : pars.momentum().unit();

  // check if you have to stop at the endSurface
  double pathLimit = std::numeric_limits<double>::infinity();
  if (endSurface) {
    // intersect the end surface
    // - it is the final one don't use the bounday check at all
    Intersection endInter
        = endSurface->intersectionEstimate(pos, dir, pDir, false);
    // non-valid intersection with the end surface provided at this layer
    // indicates wrong direction or faulty setup
    // -> do not return compatible surfaces since they may lead you on a wrong
    // navigation path
    if (endInter.valid && endInter.pathLength > 0.)
      pathLimit = endInter.pathLength;
    else
      return cSurfaces;
  }

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  if (m_approachDescriptor
      && (resolvePassive || (resolveMaterial && m_ssApproachSurfaces > 1))) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces
        = m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& aSurface : approachSurfaces) {
      if (aSurface == startSurface || aSurface == endSurface) continue;
      // we fill passive always, rest is only for material
      if (resolvePassive || aSurface->associatedMaterial())
        testCompatibleSurface(
            cSurfaces, *aSurface, pos, dir, pDir, tCheck, pathLimit);
    }
  }

  // (B) sensitive surface section
  //
  bool resolvePS = resolvePassive || resolveSensitive;
  // we have to search for if m_surfaceArray exists && either :
  // - resolvePassive is set true : records everything
  // - resolveSensitive is set true : direct request
  // - resolveMaterial is set true and sensitive structure >1

  auto crit
      = [&resolvePS, &startSurface, &endSurface](const Surface* srf) -> bool {
    bool doCollect         = resolvePS || srf->associatedMaterial();
    bool startOrEndSurface = srf == startSurface || srf == endSurface;

    return doCollect && !startOrEndSurface;
  };

  if (m_surfaceArray
      && (resolvePS || (resolveMaterial && m_ssSensitiveSurfaces > 1))) {
    // compatible test surfaces
    std::vector<const Surface*> ctestSurfaces;
    if (searchType <= 0) {

      const std::vector<const Surface*>& allTestSurfaces
          = m_surfaceArray->surfaces();
      ctestSurfaces.reserve(allTestSurfaces.size());
      std::copy_if(allTestSurfaces.begin(),
                   allTestSurfaces.end(),
                   std::back_inserter(ctestSurfaces),
                   crit);

    } else {

      const std::vector<const Surface*>& neighbors
          = m_surfaceArray->neighbors(pos);

      // Returned surfaces from SurfaceArray::neighbors() also includes
      // the nominal surface. So we can just copy here
      std::copy_if(neighbors.begin(),
                   neighbors.end(),
                   std::back_inserter(ctestSurfaces),
                   crit);
    }
    // sensitive surfaces and test them
    for (auto& ctSurface : ctestSurfaces)
      testCompatibleSurface(
          cSurfaces, *ctSurface, pos, dir, pDir, tCheck, pathLimit);

  }  // end of sensitive surfaces to exist

  // (C) this is the representing surface
  //
  // the layer surface itself is a testSurface
  const Surface* layerSurface = &surfaceRepresentation();
  // veto if it is the surface of the track parameter already
  if (layerSurface != startSurface && layerSurface != endSurface
      && (resolvePassive
          || (resolveMaterial && layerSurface->associatedMaterial()))) {
    testCompatibleSurface(
        cSurfaces, *layerSurface, pos, dir, pDir, tCheck, pathLimit);
  }

  // only sort it if the intersection was done
  std::sort(cSurfaces.begin(), cSurfaces.end());

  // return
  return cSurfaces;
}

inline void
Layer::testCompatibleSurface(std::vector<SurfaceIntersection>& cSurfaces,
                             const Surface&                    surface,
                             const Vector3D&                   pos,
                             const Vector3D&                   dir,
                             NavigationDirection               navDir,
                             const BoundaryCheck&              bcheck,
                             double                            pathLimit) const
{
  // the intersection
  Intersection sfIntersection
      = surface.intersectionEstimate(pos, dir, navDir, bcheck);
  // check if intersection is valid and pathLimit has not been exceeded
  if (sfIntersection.valid
      && sfIntersection.pathLength < pathLimit) {  // and the surfaces &
                                                   // direction to push back
                                                   // - take all
    cSurfaces.push_back(SurfaceIntersection(sfIntersection, &surface, navDir));
  }
}

inline const SurfaceIntersection
Layer::surfaceOnApproach(const Vector3D&      position,
                         const Vector3D&      momentum,
                         NavigationDirection  nDir,
                         const BoundaryCheck& bcheck,
                         bool                 resolveSensitive,
                         bool                 resolveMaterial,
                         bool                 resolvePassive) const
{
  // we need the approach surface when

  // - resolvePassive is on -> always
  // - resolveSensitive is on -> always
  // - resolveMaterial is on

  // Internal direction
  // - is nDir, but forward if anyDirection was chosen
  NavigationDirection iDir = (nDir == anyDirection) ? forward : nDir;
  //   && either sensitive or approach surfaces have material
  bool resolvePS = resolveSensitive || resolvePassive;
  bool resolveMS = resolveMaterial
      && (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1
          || surfaceRepresentation().associatedMaterial());
  // now of course this only counts when you have an approach descriptor
  if (m_approachDescriptor && (resolvePS || resolveMS)) {
    // test if you are on an approach surface already, if so - provide
    if (nDir == anyDirection) {
      for (auto& asf : m_approachDescriptor->containedSurfaces()) {
        if (asf->isOnSurface(position, bcheck)) {
          Intersection nIntersection(position, 0., true);
          return SurfaceIntersection(
              nIntersection, &surfaceRepresentation(), nDir);
        }
      }
    }
    // that's the collect trigger for always collecting
    // let's find the approach surface
    SurfaceIntersection aSurface = m_approachDescriptor->approachSurface(
        position, momentum, iDir, bcheck);
    if (aSurface.intersection.valid) return (aSurface);
  }
  // allow to stay if you are on the surface
  if (nDir == anyDirection
      && surfaceRepresentation().isOnSurface(position, bcheck)) {
    // create a valid 0-distance intescetion and return it
    Intersection nIntersection(position, 0., true);
    return SurfaceIntersection(nIntersection, &surfaceRepresentation(), nDir);
  }

  // create the intersection with the surface representation
  auto rIntersection = surfaceRepresentation().intersectionEstimate(
      position, momentum, iDir, bcheck);
  // return the result
  return SurfaceIntersection(rIntersection, &surfaceRepresentation(), iDir);
}

inline bool
Layer::isOnLayer(const Vector3D& gp, const BoundaryCheck& bcheck) const
{
  if (m_representingVolume) return m_representingVolume->inside(gp);
  return (surfaceRepresentation()).isOnSurface(gp, bcheck);
}

inline std::vector<SurfaceIntersection>
Layer::compatibleSurfaces(const TrackParameters& pars,
                          NavigationDirection    pdir,
                          const BoundaryCheck&   bcheck,
                          bool                   resolveSensitive,
                          bool                   resolveMaterial,
                          bool                   resolvePassive,
                          int                    searchType,
                          const Surface*         startSurface,
                          const Surface*         endSurface) const
{
  return getCompatibleSurfaces(pars,
                               pdir,
                               bcheck,
                               resolveSensitive,
                               resolveMaterial,
                               resolvePassive,
                               searchType,
                               startSurface,
                               endSurface);
}

inline std::vector<SurfaceIntersection>
Layer::compatibleSurfaces(const NeutralParameters& pars,
                          NavigationDirection      pdir,
                          const BoundaryCheck&     bcheck,
                          bool                     resolveSensitive,
                          bool                     resolveMaterial,
                          bool                     resolvePassive,
                          int                      searchType,
                          const Surface*           startSurface,
                          const Surface*           endSurface) const
{
  return getCompatibleSurfaces(pars,
                               pdir,
                               bcheck,
                               resolveSensitive,
                               resolveMaterial,
                               resolvePassive,
                               searchType,
                               startSurface,
                               endSurface);
}

// ----------- legacy method block: end ----------------------

}  // namespace Acts
