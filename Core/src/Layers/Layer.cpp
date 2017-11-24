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

Acts::Layer::Layer(std::unique_ptr<SurfaceArray_old>   surfaceArray,
                   double                              thickness,
                   std::unique_ptr<ApproachDescriptor> ades,
                   LayerType                           laytyp)
  : m_nextLayers(NextLayers(nullptr, nullptr))
  , m_nextLayerUtility(nullptr)
  , m_surfaceArray(surfaceArray.release())
  , m_layerThickness(thickness)
  , m_approachDescriptor(nullptr)
  , m_enclosingTrackingVolume(nullptr)
  , m_enclosingDetachedTrackingVolume(nullptr)
  , m_representingVolume(nullptr)
  , m_layerType(laytyp)
  , m_ssRepresentingSurface(1)
  , m_ssSensitiveSurfaces(0)
  , m_ssApproachSurfaces(0)
  , m_detectorElements()

{
  if (ades) {
    ades->registerLayer(*this);
    m_approachDescriptor = std::move(ades);
    m_ssApproachSurfaces = 1;  // indicates existence
  }
  // indicates existence of sensitive surfaces
  if (m_surfaceArray) m_ssSensitiveSurfaces = 1;
}

Acts::Layer::~Layer()
{
  delete m_representingVolume;
}

bool
Acts::Layer::isOnLayer(const Acts::Vector3D& gp,
                       const BoundaryCheck&  bcheck) const
{
  return (surfaceRepresentation()).isOnSurface(gp, bcheck);
}

const Acts::SurfaceIntersection
Acts::Layer::surfaceOnApproach(const Acts::Vector3D&      position,
                               const Acts::Vector3D&      momentum,
                               Acts::PropDirection        pDir,
                               const Acts::BoundaryCheck& bcheck,
                               bool                       collectSensitive,
                               bool                       collectMaterial,
                               bool                       collectPassive,
                               const Acts::ICompatibilityEstimator*) const
{
  // we need the approach surface when
  // - collectPassive is on -> always
  // - collectSensitive is on -> always
  // - collectMaterial is on
  //   && either sensitive or approach surfaces have material
  bool collectPS = collectSensitive || collectPassive;
  bool collectMS = collectMaterial
      && (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1
          || surfaceRepresentation().associatedMaterial());
  // now of course this only counts when you have an approach descriptor
  if (m_approachDescriptor && (collectPS || collectMS)) {
    // that's the collect trigger for always collecting
    // let's find the approach surface
    SurfaceIntersection aSurface = m_approachDescriptor->approachSurface(
        position, pDir * momentum.unit(), bcheck);
    if (aSurface.intersection.valid) return (aSurface);
  }
  // create the intersection with the surface reprensentation
  auto rIntersection = surfaceRepresentation().intersectionEstimate(
      position, pDir * momentum.unit(), true, bcheck);
  // return the result
  return SurfaceIntersection(rIntersection, &surfaceRepresentation(), pDir);
}

std::vector<Acts::SurfaceIntersection>
Acts::Layer::compatibleSurfaces(const Acts::TrackParameters& pars,
                                Acts::PropDirection          pdir,
                                const Acts::BoundaryCheck&   bcheck,
                                bool                         collectSensitive,
                                bool                         collectMaterial,
                                bool                         collectPassive,
                                int                          searchType,
                                const Acts::Surface*         startSurface,
                                const Acts::Surface*         endSurface,
                                const Acts::ICompatibilityEstimator* ice) const
{
  return getCompatibleSurfaces(pars,
                               pdir,
                               bcheck,
                               collectSensitive,
                               collectMaterial,
                               collectPassive,
                               searchType,
                               startSurface,
                               endSurface,
                               ice);
}

std::vector<Acts::SurfaceIntersection>
Acts::Layer::compatibleSurfaces(const Acts::NeutralParameters& pars,
                                Acts::PropDirection            pdir,
                                const Acts::BoundaryCheck&     bcheck,
                                bool                           collectSensitive,
                                bool                           collectMaterial,
                                bool                           collectPassive,
                                int                            searchType,
                                const Acts::Surface*           startSurface,
                                const Acts::Surface*           endSurface,
                                const Acts::ICompatibilityEstimator* ice) const
{
  return getCompatibleSurfaces(pars,
                               pdir,
                               bcheck,
                               collectSensitive,
                               collectMaterial,
                               collectPassive,
                               searchType,
                               startSurface,
                               endSurface,
                               ice);
}

const Acts::ApproachDescriptor*
Acts::Layer::approachDescriptor() const
{
  return m_approachDescriptor.get();
}

Acts::ApproachDescriptor*
Acts::Layer::approachDescriptor()
{
  return const_cast<ApproachDescriptor*>(m_approachDescriptor.get());
}

void
Acts::Layer::closeGeometry(const GeometryID& layerID)
{
  // set the volumeID of this
  assignGeoID(layerID);

  // also find out how the sub structure is defined
  if (surfaceRepresentation().associatedMaterial()) m_ssRepresentingSurface = 2;

  // loop over the approach surfaces
  if (m_approachDescriptor) {
    // indicates the existance of approach surfaces
    m_ssApproachSurfaces = 1;
    // loop through the approachSurfaces and assign unique GeomeryID
    geo_id_value iasurface = 0;
    for (auto& aSurface : m_approachDescriptor->containedSurfaces()) {
      GeometryID asurfaceID = layerID;
      asurfaceID.add(++iasurface, GeometryID::approach_mask);
      auto mutableASurface = const_cast<Surface*>(aSurface);
      mutableASurface->assignGeoID(asurfaceID);
      // if any of the approach surfaces has material
      if (aSurface->associatedMaterial()) m_ssApproachSurfaces = 2;
    }
  }
  // check if you have sensitive surfaces
  if (m_surfaceArray) {
    // indicates the existance of sensitive surfaces
    m_ssSensitiveSurfaces = 1;
    // loop sensitive surfaces and assign unique GeometryID
    geo_id_value issurface = 0;
    for (auto& sSurface : m_surfaceArray->arrayObjects()) {
      GeometryID ssurfaceID = layerID;
      ssurfaceID.add(++issurface, GeometryID::sensitive_mask);
      auto mutableSSurface = const_cast<Surface*>(sSurface);
      mutableSSurface->assignGeoID(ssurfaceID);
      // fill the map of detector elements
      m_detectorElements.emplace(sSurface->associatedIdentifier(),
                                 sSurface->associatedDetectorElement());
      // if any of the sensitive surfaces has material
      if (sSurface->associatedMaterial()) m_ssSensitiveSurfaces = 2;
    }
  }
}
