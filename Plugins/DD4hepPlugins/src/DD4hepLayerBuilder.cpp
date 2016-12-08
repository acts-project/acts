// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepLayerBuilder.hpp"
#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"
#include "ACTS/Plugins/DD4hepPlugins/IActsExtension.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Tools/ILayerCreator.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "DD4hep/Detector.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::DD4hepLayerBuilder::DD4hepLayerBuilder(
    const Acts::DD4hepLayerBuilder::Config& config,
    std::unique_ptr<Logger>                 logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(config);
}

Acts::DD4hepLayerBuilder::~DD4hepLayerBuilder()
{
}

void
Acts::DD4hepLayerBuilder::setConfiguration(
    const Acts::DD4hepLayerBuilder::Config& config)
{
  m_cfg = config;
}

const Acts::LayerVector
Acts::DD4hepLayerBuilder::negativeLayers() const
{
  LayerVector layers;
  if (m_cfg.negativeLayers.empty()) {
    ACTS_VERBOSE("[L] No layers handed over for negative volume.");
  } else {
    ACTS_VERBOSE("[L] Received layers for negative volume -> creating "
                 "disc layers");
    // go through layers
    for (auto& detElement : m_cfg.negativeLayers) {
      // prepare the layer surfaces
      std::vector<const Surface*> layerSurfaces;
      // access the extension of the layer
      // at this stage all layer detElements have extension (checked in
      // ConvertDD4hepDetector)
      Acts::IActsExtension* detExtension
          = detElement.extension<Acts::IActsExtension>();
      // access the axis orienation of the modules
      std::string axes = detExtension->axes();
      // collect the sensitive detector elements possibly contained by the layer
      collectSensitive(detElement, layerSurfaces, axes);
      // access the global transformation matrix of the layer
      auto transform = convertTransform(&(detElement.worldTransformation()));
      // get the shape of the layer
      TGeoShape* geoShape
          = detElement.placement().ptr()->GetVolume()->GetShape();
      if (geoShape) {
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube)
          ACTS_ERROR(
              "[L] Disc layer has wrong shape - needs to be TGeoConeSeg!");
        // extract the boundaries
        double rMin = tube->GetRmin1() * units::_cm;
        double rMax = tube->GetRmax1() * units::_cm;
        double zMin
            = (transform->translation()
               - transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        double zMax
            = (transform->translation()
               + transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        if (zMin > zMax) std::swap(zMin, zMax);
        layers.push_back(m_cfg.layerCreator->discLayer(layerSurfaces,
                                                       zMin,
                                                       zMax,
                                                       rMin,
                                                       rMax,
                                                       m_cfg.bTypeR,
                                                       m_cfg.bTypePhi,
                                                       transform));
      } else if (detExtension->buildEnvelope()) {
        layers.push_back(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          detExtension->envelopeR(),
                                          detExtension->envelopeR(),
                                          detExtension->envelopeZ(),
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));
    }
  }
  return layers;
}

const Acts::LayerVector
Acts::DD4hepLayerBuilder::centralLayers() const
{
  LayerVector layers;
  if (m_cfg.centralLayers.empty()) {
    ACTS_VERBOSE("[L] No layers handed over for central volume!");
  } else {
    ACTS_VERBOSE("[L] Received layers for central volume -> creating "
                 "cylindrical layers");
    // go through layers
    for (auto& detElement : m_cfg.centralLayers) {
      // prepare the layer surfaces
      std::vector<const Surface*> layerSurfaces;
      // access the extension of the layer
      // at this stage all layer detElements have extension (checked in
      // ConvertDD4hepDetector)
      Acts::IActsExtension* detExtension
          = detElement.extension<Acts::IActsExtension>();
      // access the axis orienation of the modules
      std::string axes = detExtension->axes();
      // collect the sensitive detector elements possibly contained by the layer
      collectSensitive(detElement, layerSurfaces, axes);
      // access the global transformation matrix of the layer
      auto transform = convertTransform(&(detElement.worldTransformation()));
      // get the shape of the layer
      TGeoShape* geoShape
          = detElement.placement().ptr()->GetVolume()->GetShape();

      if (detExtension->buildEnvelope()) {
        layers.push_back(
            m_cfg.layerCreator->cylinderLayer(layerSurfaces,
                                              detExtension->envelopeR(),
                                              detExtension->envelopeZ(),
                                              m_cfg.bTypePhi,
                                              m_cfg.bTypeZ,
                                              transform));
      } else if (geoShape) {
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube)
          ACTS_ERROR(
              "[L] Cylinder layer has wrong shape - needs to be TGeoConeSeg!");
        // extract the boundaries
        double rMin  = tube->GetRmin1() * units::_cm;
        double rMax  = tube->GetRmax1() * units::_cm;
        double halfZ = tube->GetDz() * units::_cm;
        layers.push_back(m_cfg.layerCreator->cylinderLayer(layerSurfaces,
                                                           rMin,
                                                           rMax,
                                                           halfZ,
                                                           m_cfg.bTypePhi,
                                                           m_cfg.bTypeZ,
                                                           transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));
    }
  }
  return layers;
}

const Acts::LayerVector
Acts::DD4hepLayerBuilder::positiveLayers() const
{
  LayerVector layers;
  if (m_cfg.positiveLayers.empty()) {
    ACTS_VERBOSE("[L] No layers handed over for negative volume!");
  } else {
    ACTS_VERBOSE("[L] Received layers for negative volume -> creating "
                 "disc layers");
    // go through layers
    for (auto& detElement : m_cfg.positiveLayers) {
      // prepare the layer surfaces
      std::vector<const Surface*> layerSurfaces;
      // access the extension of the layer
      // at this stage all layer detElements have extension (checked in
      // ConvertDD4hepDetector)
      Acts::IActsExtension* detExtension
          = detElement.extension<Acts::IActsExtension>();
      // access the axis orienation of the modules
      std::string axes = detExtension->axes();
      // collect the sensitive detector elements possibly contained by the layer
      collectSensitive(detElement, layerSurfaces, axes);
      // access the global transformation matrix of the layer
      auto transform = convertTransform(&(detElement.worldTransformation()));
      // get the shape of the layer
      TGeoShape* geoShape
          = detElement.placement().ptr()->GetVolume()->GetShape();
      if (geoShape) {
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube)
          ACTS_ERROR(
              "[L] Disc layer has wrong shape - needs to be TGeoConeSeg!");
        // extract the boundaries
        double rMin = tube->GetRmin1() * units::_cm;
        double rMax = tube->GetRmax1() * units::_cm;
        double zMin
            = (transform->translation()
               - transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        double zMax
            = (transform->translation()
               + transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        if (zMin > zMax) std::swap(zMin, zMax);
        layers.push_back(m_cfg.layerCreator->discLayer(layerSurfaces,
                                                       zMin,
                                                       zMax,
                                                       rMin,
                                                       rMax,
                                                       m_cfg.bTypeR,
                                                       m_cfg.bTypePhi,
                                                       transform));
      } else if (detExtension->buildEnvelope()) {
        layers.push_back(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          detExtension->envelopeR(),
                                          detExtension->envelopeR(),
                                          detExtension->envelopeZ(),
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));
    }
  }
  return layers;
}

void
Acts::DD4hepLayerBuilder::collectSensitive(
    const DD4hep::Geometry::DetElement& detElement,
    std::vector<const Acts::Surface*>&  surfaces,
    const std::string&                  axes) const
{
  const DD4hep::Geometry::DetElement::Children& children
      = detElement.children();
  if (!children.empty()) {
    for (auto& child : children) {
      DD4hep::Geometry::DetElement childDetElement = child.second;
      if (childDetElement.volume().isSensitive()) {
        // create the corresponding detector element
        Acts::DD4hepDetElement* dd4hepDetElement
            = new Acts::DD4hepDetElement(childDetElement, axes, units::_cm);
        // add surface to surface vector
        surfaces.push_back(&(dd4hepDetElement->surface()));
      }
      collectSensitive(childDetElement, surfaces);
    }
  }
}

std::shared_ptr<Acts::Transform3D>
Acts::DD4hepLayerBuilder::convertTransform(const TGeoMatrix* tGeoTrans) const
{
  // get the placement and orientation in respect to its mother
  const Double_t* rotation    = tGeoTrans->GetRotationMatrix();
  const Double_t* translation = tGeoTrans->GetTranslation();
  auto            transform   = std::make_shared<Acts::Transform3D>(
      Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
      Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
      Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
      Acts::Vector3D(translation[0] * units::_cm,
                     translation[1] * units::_cm,
                     translation[2] * units::_cm));
  return (transform);
}
