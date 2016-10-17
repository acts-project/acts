// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetachedTrackingVolume.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Detector/DetachedTrackingVolume.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"

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
  , m_trkVolume(volume)
  , m_layerRepresentation(layer)
  , m_multilayerRepresentation(multiLayer)
  , m_baseTransform(nullptr)
  , m_constituents(nullptr)
{
}

Acts::DetachedTrackingVolume::~DetachedTrackingVolume()
{
  delete m_baseTransform;
}

void
Acts::DetachedTrackingVolume::move(Acts::Transform3D&) const
{
  //!< @todo implement
  //// move the volume
  // m_trkVolume->moveTrackingVolume(shift);
  //// move the single layer
  // if (m_layerRepresentation)
  //    m_layerRepresentation->moveLayer(shift);
  //// move the multi layer
  // if (m_multilayerRepresentation)
  //    for (auto& layer: (*m_multilayerRepresentation))
  //        layer->moveLayer(shift);
}

std::shared_ptr<const Acts::DetachedTrackingVolume>
Acts::DetachedTrackingVolume::clone(std::string        name,
                                    Acts::Transform3D& shift) const
{
  // create the new base tracking volume
  std::shared_ptr<const TrackingVolume> shiftedTrackingVolume
    = TrackingVolume::create(*m_trkVolume,shift);

  // create and shift the layers if there are any
  std::shared_ptr<const Layer> layerRepresentation = m_layerRepresentation
      ? m_layerRepresentation->cloneWithShift(shift)
      : nullptr;
  std::vector<std::shared_ptr<const Layer>> multiLayerRepresentation;

  // create and shift the multi layer
  if (!m_multilayerRepresentation.empty()) {
    // create the new multi layer representation
    multiLayerRepresentation.reserve(m_multilayerRepresentation.size());
    // now clone the original consitutents
    for (auto& lay : m_multilayerRepresentation)
      multiLayerRepresentation.push_back(lay->cloneWithShift(shift));
  }

  // return a new cloned TrackingVolume
  return std::shared_ptr<const Acts::DetachedTrackingVolume>(
      new Acts::DetachedTrackingVolume(name,
                                       shiftedTrackingVolume,
                                       layerRepresentation,
                                       multiLayerRepresentation));
}

void
Acts::DetachedTrackingVolume::sign(GeometrySignature signat,
                                   GeometryType      geotype) const
{
  m_trkVolume->sign(signat, geotype);
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
Acts::DetachedTrackingVolume::setBaseTransform(Acts::Transform3D* transf) const
{
  if (transf)
    m_baseTransform = transf;
  else {
    delete m_baseTransform;
    m_baseTransform
        = new Acts::Transform3D(this->trackingVolume()->transform());
  }
}

void
Acts::DetachedTrackingVolume::realign(Acts::Transform3D* transf) const
{
  if (transf) {
    Acts::Transform3D shift
        = (*transf) * this->trackingVolume()->transform().inverse();
    this->move(shift);
  } else if (m_baseTransform) {
    Acts::Transform3D shift
        = *m_baseTransform * this->trackingVolume()->transform().inverse();
    this->move(shift);
  }
}
