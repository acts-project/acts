// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Layer.cpp, Acts project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "Acts/Layers/Layer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ApproachDescriptor.hpp"
#include "Acts/Utilities/BinUtility.hpp"

Acts::Layer::Layer(std::unique_ptr<SurfaceArray>       surfaceArray,
                   double                              thickness,
                   std::unique_ptr<ApproachDescriptor> ades,
                   LayerType                           laytyp)
  : m_nextLayers(NextLayers(nullptr, nullptr))
  , m_nextLayerUtility(nullptr)
  , m_surfaceArray(surfaceArray.release())
  , m_layerThickness(thickness)
  , m_approachDescriptor(nullptr)
  , m_trackingVolume(nullptr)
  , m_enclosingDetachedTrackingVolume(nullptr)
  , m_representingVolume(nullptr)
  , m_layerType(laytyp)
  , m_ssRepresentingSurface(1)
  , m_ssSensitiveSurfaces(0)
  , m_ssApproachSurfaces(0)

{
  if (ades) {
    ades->registerLayer(*this);
    m_approachDescriptor = std::move(ades);
    m_ssApproachSurfaces = 1;  // indicates existence
  }
  // indicates existence of sensitive surfaces
  if (m_surfaceArray) {
    m_ssSensitiveSurfaces = 1;
  }
}

Acts::Layer::~Layer()
{
  delete m_representingVolume;
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
  if (surfaceRepresentation().associatedMaterial() != nullptr) {
    m_ssRepresentingSurface = 2;
  }

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
      if (aSurface->associatedMaterial() != nullptr) {
        m_ssApproachSurfaces = 2;
      }
    }
  }
  // check if you have sensitive surfaces
  if (m_surfaceArray) {
    // indicates the existance of sensitive surfaces
    m_ssSensitiveSurfaces = 1;
    // loop sensitive surfaces and assign unique GeometryID
    geo_id_value issurface = 0;
    for (auto& sSurface : m_surfaceArray->surfaces()) {
      GeometryID ssurfaceID = layerID;
      ssurfaceID.add(++issurface, GeometryID::sensitive_mask);
      auto mutableSSurface = const_cast<Surface*>(sSurface);
      mutableSSurface->assignGeoID(ssurfaceID);
      // if any of the sensitive surfaces has material
      if (sSurface->associatedMaterial() != nullptr) {
        m_ssSensitiveSurfaces = 2;
      }
    }
  }
}
