// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepLayerBuilder.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Units.hpp"
#include "DD4hep/Detector.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include <boost/algorithm/string.hpp>

Acts::DD4hepLayerBuilder::DD4hepLayerBuilder(
    const Acts::DD4hepLayerBuilder::Config& config,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(config);
}

Acts::DD4hepLayerBuilder::~DD4hepLayerBuilder() = default;

void Acts::DD4hepLayerBuilder::setConfiguration(
    const Acts::DD4hepLayerBuilder::Config& config) {
  m_cfg = config;
}

const Acts::LayerVector Acts::DD4hepLayerBuilder::endcapLayers(
    const GeometryContext& gctx,
    const std::vector<dd4hep::DetElement>& dendcapLayers,
    const std::string& side) const {
  LayerVector layers;
  if (dendcapLayers.empty()) {
    ACTS_VERBOSE(" No layers handed over for " << side << " volume!");
  } else {
    ACTS_VERBOSE(" Received layers for " << side
                                         << " volume -> creating "
                                            "disc layers");
    // go through layers
    for (auto& detElement : dendcapLayers) {
      // prepare the layer surfaces
      std::vector<std::shared_ptr<const Surface>> layerSurfaces;
      // access the extension of the layer
      // at this stage all layer detElements have extension (checked in
      // ConvertDD4hepDetector)
      Acts::ActsExtension* detExtension =
          detElement.extension<Acts::ActsExtension>();
      // collect the sensitive detector elements possibly contained by the layer
      resolveSensitive(detElement, layerSurfaces);
      // access the global transformation matrix of the layer
      auto transform =
          convertTransform(&(detElement.nominal().worldTransformation()));
      // get the shape of the layer
      TGeoShape* geoShape =
          detElement.placement().ptr()->GetVolume()->GetShape();
      // create the proto layer
      ProtoLayer pl(gctx, layerSurfaces);
      if (detExtension->hasValue("r", "envelope") &&
          detExtension->hasValue("z", "envelope")) {
        // set the values of the proto layer in case enevelopes are handed over
        pl.envelope[Acts::binR] = {detExtension->getValue("r", "envelope"),
                                   detExtension->getValue("r", "envelope")};
        pl.envelope[Acts::binZ] = {detExtension->getValue("z", "envelope"),
                                   detExtension->getValue("z", "envelope")};
      } else if (geoShape != nullptr) {
        TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(geoShape);
        if (tube == nullptr) {
          ACTS_ERROR(" Disc layer has wrong shape - needs to be TGeoTubeSeg!");
        }
        // extract the boundaries
        double rMin = tube->GetRmin() * UnitConstants::cm;
        double rMax = tube->GetRmax() * UnitConstants::cm;
        double zMin =
            (transform->translation() -
             transform->rotation().col(2) * tube->GetDz() * UnitConstants::cm)
                .z();
        double zMax =
            (transform->translation() +
             transform->rotation().col(2) * tube->GetDz() * UnitConstants::cm)
                .z();
        if (zMin > zMax) {
          std::swap(zMin, zMax);
        }
        // check if layer has surfaces
        if (layerSurfaces.empty()) {
          ACTS_VERBOSE(" Disc layer has no senstive surfaces.");
          // in case no surfaces are handed over the layer thickness will be set
          // to a default value to allow attaching material layers
          double z = (zMin + zMax) * 0.5;
          // create layer without surfaces
          // manually create a proto layer
          double eiz = (z != 0.) ? z - m_cfg.defaultThickness : 0.;
          double eoz = (z != 0.) ? z + m_cfg.defaultThickness : 0.;
          pl.extent.ranges[Acts::binZ] = {eiz, eoz};
          pl.extent.ranges[Acts::binR] = {rMin, rMax};
          pl.envelope[Acts::binR] = {0., 0.};
          pl.envelope[Acts::binZ] = {0., 0.};
        } else {
          ACTS_VERBOSE(" Disc layer has " << layerSurfaces.size()
                                          << " senstive surfaces.");
          // set the values of the proto layer in case dimensions are given by
          // geometry
          pl.envelope[Acts::binZ] = {std::abs(zMin - pl.min(Acts::binZ)),
                                     std::abs(zMax - pl.max(Acts::binZ))};
          pl.envelope[Acts::binR] = {std::abs(rMin - pl.min(Acts::binR)),
                                     std::abs(rMax - pl.max(Acts::binR))};
        }
      } else {
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name() +
            std::string(" has neither a shape nor tolerances for envelopes "
                        "added to its extension. Please check your detector "
                        "constructor!"));
      }

      std::shared_ptr<Layer> endcapLayer = nullptr;
      // In case the layer is sensitive
      if (detElement.volume().isSensitive()) {
        // Create the sensitive surface
        auto sensitiveSurf = createSensitiveSurface(detElement, true);
        // Create the surfaceArray
        std::unique_ptr<Acts::SurfaceArray> sArray =
            std::make_unique<SurfaceArray>(sensitiveSurf);

        // create the share disc bounds
        auto dBounds = std::make_shared<const RadialBounds>(pl.min(Acts::binR),
                                                            pl.max(Acts::binR));
        double thickness = std::fabs(pl.max(Acts::binZ) - pl.min(Acts::binZ));
        // Create the layer containing the sensitive surface
        endcapLayer = DiscLayer::create(transform, dBounds, std::move(sArray),
                                        thickness, nullptr, Acts::active);

      } else {
        endcapLayer = m_cfg.layerCreator->discLayer(
            gctx, layerSurfaces, m_cfg.bTypeR, m_cfg.bTypePhi, pl, transform,
            nullptr);
      }
      // Add the ProtoMaterial if present
      addDiscLayerProtoMaterial(detElement, *endcapLayer);
      // push back created layer
      layers.push_back(endcapLayer);
    }
  }
  return layers;
}

const Acts::LayerVector Acts::DD4hepLayerBuilder::negativeLayers(
    const GeometryContext& gctx) const {
  return endcapLayers(gctx, m_cfg.negativeLayers, "negative");
}

const Acts::LayerVector Acts::DD4hepLayerBuilder::centralLayers(
    const GeometryContext& gctx) const {
  LayerVector layers;
  if (m_cfg.centralLayers.empty()) {
    ACTS_VERBOSE(" No layers handed over for central volume!");
  } else {
    ACTS_VERBOSE(
        " Received layers for central volume -> creating "
        "cylindrical layers");
    // go through layers
    for (auto& detElement : m_cfg.centralLayers) {
      // prepare the layer surfaces
      std::vector<std::shared_ptr<const Surface>> layerSurfaces;
      // access the extension of the layer
      // at this stage all layer detElements have extension (checked in
      // ConvertDD4hepDetector)
      Acts::ActsExtension* detExtension =
          detElement.extension<Acts::ActsExtension>();
      // collect the sensitive detector elements possibly contained by the layer
      resolveSensitive(detElement, layerSurfaces);
      // access the global transformation matrix of the layer
      auto transform =
          convertTransform(&(detElement.nominal().worldTransformation()));
      // get the shape of the layer
      TGeoShape* geoShape =
          detElement.placement().ptr()->GetVolume()->GetShape();
      // create the proto layer
      ProtoLayer pl(gctx, layerSurfaces);
      if (detExtension->hasValue("r", "envelope") &&
          detExtension->hasValue("z", "envelope")) {
        // set the values of the proto layer in case enevelopes are handed over
        pl.envelope[Acts::binR] = {detExtension->getValue("r", "envelope"),
                                   detExtension->getValue("r", "envelope")};
        pl.envelope[Acts::binZ] = {detExtension->getValue("z", "envelope"),
                                   detExtension->getValue("z", "envelope")};
      } else if (geoShape != nullptr) {
        TGeoTubeSeg* tube = dynamic_cast<TGeoTubeSeg*>(geoShape);
        if (tube == nullptr)
          ACTS_ERROR(
              " Cylinder layer has wrong shape - needs to be TGeoTubeSeg!");

        // extract the boundaries
        double rMin = tube->GetRmin() * UnitConstants::cm;
        double rMax = tube->GetRmax() * UnitConstants::cm;
        double dz = tube->GetDz() * UnitConstants::cm;
        // check if layer has surfaces
        if (layerSurfaces.empty()) {
          // in case no surfaces are handed over the layer thickness will be set
          // to a default value to allow attaching material layers
          double r = (rMin + rMax) * 0.5;
          // create layer without surfaces
          // manually create a proto layer
          double eir = (r != 0.) ? r - m_cfg.defaultThickness : 0.;
          double eor = (r != 0.) ? r + m_cfg.defaultThickness : 0.;
          pl.extent.ranges[Acts::binR] = {eir, eor};
          pl.extent.ranges[Acts::binZ] = {-dz, dz};
          pl.envelope[Acts::binR] = {0., 0.};
          pl.envelope[Acts::binZ] = {0., 0.};
        } else {
          // set the values of the proto layer in case dimensions are given by
          // geometry
          pl.envelope[Acts::binZ] = {std::abs(-dz - pl.min(Acts::binZ)),
                                     std::abs(dz - pl.max(Acts::binZ))};
          pl.envelope[Acts::binR] = {std::abs(rMin - pl.min(Acts::binR)),
                                     std::abs(rMax - pl.max(Acts::binR))};
        }
      } else {
        throw std::logic_error(
            std::string("Layer DetElement: ") + detElement.name() +
            std::string(" has neither a shape nor tolerances for envelopes "
                        "added to it¥s extension. Please check your detector "
                        "constructor!"));
      }

      double halfZ = (pl.min(Acts::binZ) - pl.max(Acts::binZ)) * 0.5;

      std::shared_ptr<Layer> centralLayer = nullptr;
      // In case the layer is sensitive
      if (detElement.volume().isSensitive()) {
        // Create the sensitive surface
        auto sensitiveSurf = createSensitiveSurface(detElement);
        // Create the surfaceArray
        std::unique_ptr<Acts::SurfaceArray> sArray =
            std::make_unique<SurfaceArray>(sensitiveSurf);

        // create the layer
        double layerR = (pl.min(Acts::binR) + pl.max(Acts::binR)) * 0.5;
        double thickness = std::fabs(pl.max(Acts::binR) - pl.min(Acts::binR));
        std::shared_ptr<const CylinderBounds> cBounds(
            new CylinderBounds(layerR, halfZ));
        // Create the layer containing the sensitive surface
        centralLayer =
            CylinderLayer::create(transform, cBounds, std::move(sArray),
                                  thickness, nullptr, Acts::active);

      } else {
        centralLayer = m_cfg.layerCreator->cylinderLayer(
            gctx, layerSurfaces, m_cfg.bTypePhi, m_cfg.bTypeZ, pl, transform,
            nullptr);
      }
      // Add the ProtoMaterial if present
      addCylinderLayerProtoMaterial(detElement, *centralLayer);
      // push back created layer
      layers.push_back(centralLayer);
    }
  }
  return layers;
}

const Acts::LayerVector Acts::DD4hepLayerBuilder::positiveLayers(
    const GeometryContext& gctx) const {
  return endcapLayers(gctx, m_cfg.positiveLayers, "positive");
}

void Acts::DD4hepLayerBuilder::resolveSensitive(
    const dd4hep::DetElement& detElement,
    std::vector<std::shared_ptr<const Acts::Surface>>& surfaces) const {
  const dd4hep::DetElement::Children& children = detElement.children();
  if (!children.empty()) {
    for (auto& child : children) {
      dd4hep::DetElement childDetElement = child.second;
      if (childDetElement.volume().isSensitive()) {
        // create the surface
        surfaces.push_back(createSensitiveSurface(childDetElement, false));
      }
      resolveSensitive(childDetElement, surfaces);
    }
  }
}

std::shared_ptr<const Acts::Surface>
Acts::DD4hepLayerBuilder::createSensitiveSurface(
    const dd4hep::DetElement& detElement, bool isDisc) const {
  // access the possible extension of the DetElement
  Acts::ActsExtension* detExtension = nullptr;
  try {
    detExtension = detElement.extension<Acts::ActsExtension>();
  } catch (std::runtime_error& e) {
    ACTS_WARNING("Could not get Acts::Extension");
    return nullptr;
  }

  auto detAxis = detExtension->getType("axes", "definitions");
  // Create the corresponding detector element !- memory leak --!
  Acts::DD4hepDetectorElement* dd4hepDetElement =
      new Acts::DD4hepDetectorElement(detElement, detAxis, UnitConstants::cm,
                                      isDisc, nullptr, nullptr);

  // return the surface
  return dd4hepDetElement->surface().getSharedPtr();
}

std::shared_ptr<const Acts::Transform3D>
Acts::DD4hepLayerBuilder::convertTransform(const TGeoMatrix* tGeoTrans) const {
  // get the placement and orientation in respect to its mother
  const Double_t* rotation = tGeoTrans->GetRotationMatrix();
  const Double_t* translation = tGeoTrans->GetTranslation();
  auto transform =
      std::make_shared<const Transform3D>(TGeoPrimitivesHelper::makeTransform(
          Acts::Vector3D(rotation[0], rotation[3], rotation[6]),
          Acts::Vector3D(rotation[1], rotation[4], rotation[7]),
          Acts::Vector3D(rotation[2], rotation[5], rotation[8]),
          Acts::Vector3D(translation[0] * UnitConstants::cm,
                         translation[1] * UnitConstants::cm,
                         translation[2] * UnitConstants::cm)));
  return (transform);
}
