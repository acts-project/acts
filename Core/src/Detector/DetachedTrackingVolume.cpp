// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetachedTrackingVolume.cpp, Acts project
///////////////////////////////////////////////////////////////////

// Geometry module
#include <utility>

#include "Acts/Detector/DetachedTrackingVolume.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/Layer.hpp"

Acts::DetachedTrackingVolume::DetachedTrackingVolume()
  : m_name("undefined")
  , m_trkVolume()
  , m_layerRepresentation(nullptr)
  , m_multilayerRepresentation()
  , m_baseTransform(nullptr)
  , m_constituents(nullptr)
{
}

Acts::DetachedTrackingVolume::DetachedTrackingVolume(
    const std::string&                              name,
    std::shared_ptr<const Acts::TrackingVolume>     volume,
    std::shared_ptr<const Acts::Layer>              layer,
    std::vector<std::shared_ptr<const Acts::Layer>> multiLayer)
  : m_name(name)
  , m_trkVolume(std::move(volume))
  , m_layerRepresentation(std::move(layer))
  , m_multilayerRepresentation(std::move(multiLayer))
  , m_baseTransform(nullptr)
  , m_constituents(nullptr)
{
}

Acts::DetachedTrackingVolume::~DetachedTrackingVolume()
{
  delete m_baseTransform;
}

void
Acts::DetachedTrackingVolume::sign(GeometrySignature signat,
                                   GeometryType      geotype)
{
  auto mutableTrkVolume = std::const_pointer_cast<TrackingVolume>(m_trkVolume);
  mutableTrkVolume->sign(signat, geotype);
}

Acts::GeometrySignature
Acts::DetachedTrackingVolume::geometrySignature() const
{
  return m_trkVolume->geometrySignature();
}

Acts::GeometryType
Acts::DetachedTrackingVolume::geometryType() const
{
  return m_trkVolume->geometryType();
}

void
Acts::DetachedTrackingVolume::setBaseTransform(const Acts::Transform3D* transf)
{
  if (transf != nullptr) {
    m_baseTransform = transf;
  } else {
    delete m_baseTransform;
    m_baseTransform
        = new Acts::Transform3D(this->trackingVolume()->transform());
  }
}
