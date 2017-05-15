// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepLayerBuilder.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Material/SurfaceMaterialProxy.hpp"
#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"
#include "ACTS/Plugins/DD4hepPlugins/IActsExtension.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Tools/ILayerCreator.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "DD4hep/Detector.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::DD4hepLayerBuilder::DD4hepLayerBuilder(
    const Acts::DD4hepLayerBuilder::Config& config,
    std::unique_ptr<const Logger>           logger)
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
      // for the boundaries of the layer
      std::shared_ptr<const DiscLayer> layer = nullptr;

      if (detExtension->buildEnvelope()) {

        layer = std::dynamic_pointer_cast<const DiscLayer>(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          detExtension->envelopeR(),
                                          detExtension->envelopeR(),
                                          detExtension->envelopeZ(),
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else if (geoShape) {
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube)
          ACTS_ERROR(
              "[L] Disc layer has wrong shape - needs to be TGeoConeSeg!");
        // extract the boundaries
        double zMin
            = (transform->translation()
               - transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        double zMax
            = (transform->translation()
               + transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        if (zMin > zMax) std::swap(zMin, zMax);
        layer = std::dynamic_pointer_cast<const DiscLayer>(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          zMin,
                                          zMax,
                                          tube->GetRmin1() * units::_cm,
                                          tube->GetRmax1() * units::_cm,
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));

      auto discBounds = std::dynamic_pointer_cast<const RadialBounds>(
          std::shared_ptr<const SurfaceBounds>(
              layer->surfaceRepresentation().bounds().clone()));

      double rMin = discBounds->rMin();
      double rMax = discBounds->rMax();
      double zMin = (transform->translation()
                     - transform->rotation().col(2) * layer->thickness() * 0.5)
                        .z();
      double zMax = (transform->translation()
                     + transform->rotation().col(2) * layer->thickness() * 0.5)
                        .z();

      // create the two dimensional BinUtility for the material map of the layer
      // if the layer should carry material it will be marked by assigning a
      // SurfaceMaterialProxy
      std::shared_ptr<const SurfaceMaterialProxy> materialProxy(nullptr);
      // the approachdescriptor telling where the material sits on the layer
      // (inner, middle, outer) Surface
      std::unique_ptr<Acts::ApproachDescriptor> approachDescriptor = nullptr;
      // material position on the layer canbe inner, outer or center and will
      // be accessed from the ActsExtensions
      Acts::LayerMaterialPos layerPos = LayerMaterialPos::inner;
      // check if layer should have material
      if (detExtension->hasSupportMaterial()) {
        std::pair<size_t, size_t> materialBins = detExtension->materialBins();
        size_t bins1    = materialBins.first;
        size_t bins2    = materialBins.second;
        Acts::BinUtility materialBinUtil(
            bins1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
        materialBinUtil += Acts::BinUtility(
            bins2, rMin, rMax, Acts::open, Acts::binR, transform);
        // and create material proxy to mark layer for material mapping
        materialProxy
            = std::make_shared<const SurfaceMaterialProxy>(materialBinUtil);
        // access the material position
        layerPos = detExtension->layerMaterialPosition();
        ACTS_VERBOSE(
            "[L] Layer is marked to carry support material on Surface ( "
            "inner=0 / center=1 / outer=2 ) :   "
            << layerPos
            << "    with binning: ["
            << bins1
            << ", "
            << bins2
            << "]");
        // Create an approachdescriptor for the layer
        // create the new surfaces for the approachdescriptor
        std::vector<const Acts::Surface*> aSurfaces;
        // create the inner and outer boundary surfaces
        // first create the positions
        Vector3D innerPos = transform->translation()
            - transform->rotation().col(2) * layer->thickness() * 0.5;
        Vector3D outerPos = transform->translation()
            + transform->rotation().col(2) * layer->thickness() * 0.5;

        if (innerPos.z() < outerPos.z()) std::swap(innerPos, outerPos);

        zMin = innerPos.z();
        zMax = outerPos.z();

        Acts::DiscSurface* innerBoundary
            = new Acts::DiscSurface(std::make_shared<const Transform3D>(
                                        transform->rotation(), innerPos),
                                    rMin,
                                    rMax);

        Acts::DiscSurface* outerBoundary
            = new Acts::DiscSurface(std::make_shared<const Transform3D>(
                                        transform->rotation(), outerPos),
                                    rMin,
                                    rMax);

        Acts::DiscSurface* centralSurface
            = new Acts::DiscSurface(transform, rMin, rMax);

        // set material surface
        if (layerPos == Acts::LayerMaterialPos::inner)
          innerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::outer)
          outerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::central)
          centralSurface->setAssociatedMaterial(materialProxy);

        // collect approach surfaces
        aSurfaces.push_back(innerBoundary);
        aSurfaces.push_back(outerBoundary);
        aSurfaces.push_back(centralSurface);
        // create an ApproachDescriptor with standard surfaces - these
        // will be deleted by the approach descriptor
        approachDescriptor
            = std::make_unique<Acts::GenericApproachDescriptor<Acts::Surface>>(
                aSurfaces);
      }

      auto negativeLayer
          = m_cfg.layerCreator->discLayer(layerSurfaces,
                                          zMin,
                                          zMax,
                                          rMin,
                                          rMax,
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform,
                                          std::move(approachDescriptor));
      // push back created layer
      layers.push_back(negativeLayer);
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
      // the boundaries of the layer

      std::shared_ptr<const CylinderLayer> layer = nullptr;

      if (detExtension->buildEnvelope()) {
        layer = std::dynamic_pointer_cast<const CylinderLayer>(
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
        layer = std::dynamic_pointer_cast<const CylinderLayer>(
            m_cfg.layerCreator->cylinderLayer(layerSurfaces,
                                              tube->GetRmin1() * units::_cm,
                                              tube->GetRmax1() * units::_cm,
                                              tube->GetDz() * units::_cm,
                                              m_cfg.bTypePhi,
                                              m_cfg.bTypeZ,
                                              transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));
      // get the boundaries needed to build the approach descriptor
      double rMin = layer->surfaceRepresentation().bounds().r()
          - layer->thickness() * 0.5;
      double rMax = layer->surfaceRepresentation().bounds().r()
          + layer->thickness() * 0.5;
      double halfZ = layer->surfaceRepresentation().bounds().halflengthZ();

      // create the two dimensional BinUtility for the material map of the layer
      Acts::BinUtility* materialBinUtil = nullptr;
      // if the layer should carry material it will be marked by assigning a
      // SurfaceMaterialProxy
      std::shared_ptr<const SurfaceMaterialProxy> materialProxy(nullptr);
      // the approachdescriptor telling where the material sits on the layer
      // (inner, middle, outer) Surface
      std::unique_ptr<Acts::ApproachDescriptor> approachDescriptor = nullptr;
      // material position on the layer can be inner, outer or center and will
      // be accessed from the ActsExtensions
      Acts::LayerMaterialPos layerPos = LayerMaterialPos::inner;
      // check if layer should have material
      if (detExtension->hasSupportMaterial()) {
        std::pair<size_t, size_t> materialBins = detExtension->materialBins();
        size_t bins1    = materialBins.first;
        size_t bins2    = materialBins.second;
        materialBinUtil = new Acts::BinUtility(
            bins1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
        (*materialBinUtil) += Acts::BinUtility(
            bins2, -halfZ, halfZ, Acts::open, Acts::binZ, transform);
        // and create material proxy to mark layer for material mapping
        materialProxy
            = std::make_shared<const SurfaceMaterialProxy>(*materialBinUtil);
        // access the material position
        layerPos = detExtension->layerMaterialPosition();
        ACTS_VERBOSE(
            "[L] Layer is marked to carry support material on Surface ( "
            "inner=0 / center=1 / outer=2 ) :   "
            << layerPos
            << "    with binning: ["
            << bins1
            << ", "
            << bins2
            << "]");
        // Create an approachdescriptor for the layer
        // create the new surfaces for the approachdescriptor
        std::vector<const Acts::Surface*> aSurfaces;
        // create the inner boundary surface
        Acts::CylinderSurface* innerBoundary
            = new Acts::CylinderSurface(transform, rMin, halfZ);
        // create outer boundary surface
        Acts::CylinderSurface* outerBoundary
            = new Acts::CylinderSurface(transform, rMax, halfZ);
        // create the central surface
        Acts::CylinderSurface* centralSurface
            = new Acts::CylinderSurface(transform, (rMin + rMax) * 0.5, halfZ);

        // check if the material should be set to the inner or outer boundary
        // and set it in case
        if (layerPos == Acts::LayerMaterialPos::inner)
          innerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::outer)
          outerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::central)
          centralSurface->setAssociatedMaterial(materialProxy);

        // collect the surfaces
        aSurfaces.push_back(innerBoundary);
        aSurfaces.push_back(centralSurface);
        aSurfaces.push_back(outerBoundary);
        // create an ApproachDescriptor with standard surfaces - these
        // will be deleted by the approach descriptor
        approachDescriptor
            = std::make_unique<Acts::GenericApproachDescriptor<Acts::Surface>>(
                aSurfaces);
      }

      auto centralLayer
          = m_cfg.layerCreator->cylinderLayer(layerSurfaces,
                                              rMin,
                                              rMax,
                                              halfZ,
                                              m_cfg.bTypePhi,
                                              m_cfg.bTypeZ,
                                              transform,
                                              std::move(approachDescriptor));

      // push back created layer
      layers.push_back(centralLayer);
    }
  }
  return layers;
}

const Acts::LayerVector
Acts::DD4hepLayerBuilder::positiveLayers() const
{
  LayerVector layers;
  if (m_cfg.positiveLayers.empty()) {
    ACTS_VERBOSE("[L] No layers handed over for positive volume!");
  } else {
    ACTS_VERBOSE("[L] Received layers for positive volume -> creating "
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
      // for the boundaries of the layer
      std::shared_ptr<const DiscLayer> layer = nullptr;

      if (detExtension->buildEnvelope()) {
        layer = std::dynamic_pointer_cast<const DiscLayer>(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          detExtension->envelopeR(),
                                          detExtension->envelopeR(),
                                          detExtension->envelopeZ(),
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else if (geoShape) {
        TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
        if (!tube)
          ACTS_ERROR(
              "[L] Disc layer has wrong shape - needs to be TGeoConeSeg!");
        // extract the boundaries
        double zMin
            = (transform->translation()
               - transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        double zMax
            = (transform->translation()
               + transform->rotation().col(2) * tube->GetDz() * units::_cm)
                  .z();
        if (zMin > zMax) std::swap(zMin, zMax);

        layer = std::dynamic_pointer_cast<const DiscLayer>(
            m_cfg.layerCreator->discLayer(layerSurfaces,
                                          zMin,
                                          zMax,
                                          tube->GetRmin1() * units::_cm,
                                          tube->GetRmax1() * units::_cm,
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform));
      } else
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name()
            + std::string(" has neither a shape nor tolerances for envelopes "
                          "added to it¥s extension. Please check your detector "
                          "constructor!"));

      auto discBounds = std::dynamic_pointer_cast<const RadialBounds>(
          std::shared_ptr<const SurfaceBounds>(
              layer->surfaceRepresentation().bounds().clone()));

      double rMin = discBounds->rMin();
      double rMax = discBounds->rMax();
      double zMin = (transform->translation()
                     - transform->rotation().col(2) * layer->thickness() * 0.5)
                        .z();
      double zMax = (transform->translation()
                     + transform->rotation().col(2) * layer->thickness() * 0.5)
                        .z();

      // create the two dimensional BinUtility for the material map of the layer
      Acts::BinUtility* materialBinUtil = nullptr;
      // if the layer should carry material it will be marked by assigning a
      // SurfaceMaterialProxy
      std::shared_ptr<const SurfaceMaterialProxy> materialProxy(nullptr);
      // the approachdescriptor telling where the material sits on the layer
      // (inner, middle, outer) Surface
      std::unique_ptr<Acts::ApproachDescriptor> approachDescriptor = nullptr;
      // material position on the layer can be inner, outer or center and will
      // be accessed from the ActsExtensions
      Acts::LayerMaterialPos layerPos = LayerMaterialPos::inner;
      // check if layer should have material
      if (detExtension->hasSupportMaterial()) {
        std::pair<size_t, size_t> materialBins = detExtension->materialBins();
        size_t bins1    = materialBins.first;
        size_t bins2    = materialBins.second;
        materialBinUtil = new Acts::BinUtility(
            bins1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
        (*materialBinUtil) += Acts::BinUtility(
            bins2, rMin, rMax, Acts::open, Acts::binR, transform);
        // and create material proxy to mark layer for material mapping
        materialProxy
            = std::make_shared<const SurfaceMaterialProxy>(*materialBinUtil);
        // access the material position
        layerPos = detExtension->layerMaterialPosition();
        ACTS_VERBOSE(
            "[L] Layer is marked to carry support material on Surface ( "
            "inner=0 / center=1 / outer=2 ) :   "
            << layerPos
            << "    with binning: ["
            << bins1
            << ", "
            << bins2
            << "]");
        // Create an approachdescriptor for the layer
        // create the new surfaces for the approachdescriptor
        std::vector<const Acts::Surface*> aSurfaces;
        // create the inner and outer boundary surfaces
        // first create the positions
        Vector3D innerPos = transform->translation()
            - transform->rotation().col(2) * layer->thickness() * 0.5;
        Vector3D outerPos = transform->translation()
            + transform->rotation().col(2) * layer->thickness() * 0.5;

        if (innerPos.z() > outerPos.z()) std::swap(innerPos, outerPos);

        zMin = innerPos.z();
        zMax = outerPos.z();

        Acts::DiscSurface* innerBoundary
            = new Acts::DiscSurface(std::make_shared<const Transform3D>(
                                        transform->rotation(), innerPos),
                                    rMin,
                                    rMax);

        Acts::DiscSurface* outerBoundary
            = new Acts::DiscSurface(std::make_shared<const Transform3D>(
                                        transform->rotation(), outerPos),
                                    rMin,
                                    rMax);

        Acts::DiscSurface* centralSurface
            = new Acts::DiscSurface(transform, rMin, rMax);

        // set material surface
        if (layerPos == Acts::LayerMaterialPos::inner)
          innerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::outer)
          outerBoundary->setAssociatedMaterial(materialProxy);

        if (layerPos == Acts::LayerMaterialPos::central)
          centralSurface->setAssociatedMaterial(materialProxy);
        // collect approach surfaces
        aSurfaces.push_back(innerBoundary);
        aSurfaces.push_back(centralSurface);
        aSurfaces.push_back(outerBoundary);
        // create an ApproachDescriptor with standard surfaces - these
        // will be deleted by the approach descriptor
        approachDescriptor
            = std::make_unique<Acts::GenericApproachDescriptor<Acts::Surface>>(
                aSurfaces);
      }

      auto positiveLayer
          = m_cfg.layerCreator->discLayer(layerSurfaces,
                                          zMin,
                                          zMax,
                                          rMin,
                                          rMax,
                                          m_cfg.bTypeR,
                                          m_cfg.bTypePhi,
                                          transform,
                                          std::move(approachDescriptor));
      // push back created layer
      layers.push_back(positiveLayer);
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
      collectSensitive(childDetElement, surfaces, axes);
    }
  }
}

std::shared_ptr<const Acts::Transform3D>
Acts::DD4hepLayerBuilder::convertTransform(const TGeoMatrix* tGeoTrans) const
{
  // get the placement and orientation in respect to its mother
  const Double_t* rotation    = tGeoTrans->GetRotationMatrix();
  const Double_t* translation = tGeoTrans->GetTranslation();
  auto            transform   = std::make_shared<const Transform3D>(
      Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
      Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
      Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
      Acts::Vector3D(translation[0] * units::_cm,
                     translation[1] * units::_cm,
                     translation[2] * units::_cm));
  return (transform);
}
