// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Layer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/ApproachDescriptor.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

Acts::Layer::Layer()
  : m_nextLayers(NextLayers(nullptr, nullptr))
  , m_nextLayerUtility(nullptr)
  , m_surfaceArray(nullptr)
  , m_layerThickness(0.)
  , m_overlapDescriptor(nullptr)
  , m_approachDescriptor(nullptr)
  , m_enclosingTrackingVolume(nullptr)
  , m_enclosingDetachedTrackingVolume(nullptr)
  , m_representingVolume(nullptr)
  , m_layerType(0)
{
  // @TODO temporary - until GeoID service is in place
  assignGeoID(GeometryID(1));
}

Acts::Layer::Layer(std::unique_ptr<Acts::SurfaceArray> surfaceArray,
                   double                              thickness,
                   Acts::OverlapDescriptor*            olap,
                   Acts::ApproachDescriptor*           ades,
                   int                                 laytyp)
  : m_nextLayers(NextLayers(nullptr, nullptr))
  , m_nextLayerUtility(nullptr)
  , m_surfaceArray(std::move(surfaceArray))
  , m_layerThickness(thickness)
  , m_overlapDescriptor(olap)
  , m_approachDescriptor(ades)
  , m_enclosingTrackingVolume(nullptr)
  , m_enclosingDetachedTrackingVolume(nullptr)
  , m_representingVolume(nullptr)
  , m_layerType(laytyp)
{
  // @TODO temporary - until GeoID service is in place
  assignGeoID(GeometryID(1));
}

Acts::Layer::Layer(const Acts::Layer& lay)
  : m_nextLayers(NextLayers(nullptr, nullptr))
  , m_nextLayerUtility(nullptr)
  , m_surfaceArray()
  , m_layerThickness(lay.m_layerThickness)
  , m_overlapDescriptor(nullptr)
  , m_approachDescriptor(nullptr)
  , m_enclosingTrackingVolume(nullptr)
  , m_enclosingDetachedTrackingVolume(nullptr)
  , m_representingVolume(lay.m_representingVolume)
  , m_layerType(lay.m_layerType)
{
}

Acts::Layer::~Layer()
{
  delete m_overlapDescriptor;
  delete m_approachDescriptor;
  delete m_representingVolume;
}

bool
Acts::Layer::isOnLayer(const Acts::Vector3D& gp,
                       const BoundaryCheck&  bchk) const
{
  return (surfaceRepresentation()).isOnSurface(gp, bchk);
}

/** Surface seen on approach - if not defined differently, it is the
 * surfaceRepresentation() */
const Acts::SurfaceIntersection
Acts::Layer::surfaceOnApproach(const Acts::Vector3D&      position,
                               const Acts::Vector3D&      momentum,
                               Acts::PropDirection        pDir,
                               const Acts::BoundaryCheck& bcheck,
                               bool,
                               const Acts::ICompatibilityEstimator*) const
{
  // forward what the approach descriptor can give you
  if (m_approachDescriptor) {
    SurfaceIntersection aSurface = m_approachDescriptor->approachSurface(
        position, pDir * momentum.unit(), bcheck);
    if (aSurface.intersection.valid) return (aSurface);
  }
  // otherwise just return the surfaceRepresentation
  return SurfaceIntersection(Intersection(), &surfaceRepresentation(), pDir);
}

bool
Acts::Layer::compatibleSurfaces(
    std::vector<Acts::SurfaceIntersection>& cSurfaces,
    const Acts::TrackParameters&            pars,
    Acts::PropDirection                     pdir,
    const Acts::BoundaryCheck&              bcheck,
    bool                                    collectSensitive,
    bool                                    collectPassive,
    int                                     searchType,
    const Acts::Surface*                    startSurface,
    const Acts::Surface*                    endSurface,
    const Acts::ICompatibilityEstimator*    ice) const
{
  return getCompatibleSurfaces(cSurfaces,
                               pars,
                               pdir,
                               bcheck,
                               collectSensitive,
                               collectPassive,
                               searchType,
                               startSurface,
                               endSurface,
                               ice);
}

bool
Acts::Layer::compatibleSurfaces(
    std::vector<Acts::SurfaceIntersection>& cSurfaces,
    const Acts::NeutralParameters&          pars,
    Acts::PropDirection                     pdir,
    const Acts::BoundaryCheck&              bcheck,
    bool                                    collectSensitive,
    bool                                    collectPassive,
    int                                     searchType,
    const Acts::Surface*                    startSurface,
    const Acts::Surface*                    endSurface,
    const Acts::ICompatibilityEstimator*    ice) const
{
  return getCompatibleSurfaces(cSurfaces,
                               pars,
                               pdir,
                               bcheck,
                               collectSensitive,
                               collectPassive,
                               searchType,
                               startSurface,
                               endSurface,
                               ice);
}

bool
Acts::Layer::hasSubStructure(bool resolveSensitive) const
{
  if (resolveSensitive && m_surfaceArray) return true;
  return false;
}

bool
Acts::Layer::hasSensitive() const
{
  return bool(m_surfaceArray);
}

bool
Acts::Layer::hasMaterial() const
{
  return true;
}

const Acts::ApproachDescriptor*
Acts::Layer::approachDescriptor() const
{
  return m_approachDescriptor;
}
