// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "../DetUtils.h"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using dd4hep::DetElement;
using dd4hep::PlacedVolume;
using dd4hep::Volume;
using dd4hep::xml::Dimension;

namespace det {
static dd4hep::Ref_t createTkLayoutTrackerBarrel(
    dd4hep::Detector& lcdd, dd4hep::xml::Handle_t xmlElement,
    dd4hep::SensitiveDetector sensDet) {
  // shorthands
  dd4hep::xml::DetElement xmlDet =
      static_cast<dd4hep::xml::DetElement>(xmlElement);
  Dimension dimensions(xmlDet.dimensions());
  // get sensitive detector type from xml
  dd4hep::xml::Dimension sdTyp = xmlElement.child(_Unicode(sensitive));
  // sensitive detector used for all sensitive parts of this detector
  sensDet.setType(sdTyp.typeStr());

  // definition of top volume
  // has min/max dimensions of tracker for visualization etc.
  std::string detectorName = xmlDet.nameStr();
  DetElement topDetElement(detectorName, xmlDet.id());
  // detElement owns extension
  Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
  detWorldExt->addType("barrel", "detector");
  topDetElement.addExtension<Acts::ActsExtension>(detWorldExt);
  dd4hep::Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(),
                              (dimensions.zmax() - dimensions.zmin()) * 0.5);
  Volume topVolume(detectorName, topVolumeShape, lcdd.air());
  topVolume.setVisAttributes(lcdd.invisible());

  // counts all layers - incremented in the inner loop over repeat - tags
  unsigned int layerCounter = 0;
  double integratedModuleComponentThickness = 0;
  double phi = 0;
  // loop over 'layer' nodes in xml
  dd4hep::xml::Component xLayers = xmlElement.child(_Unicode(layers));
  for (dd4hep::xml::Collection_t xLayerColl(xLayers, _U(layer));
       nullptr != xLayerColl; ++xLayerColl) {
    dd4hep::xml::Component xLayer =
        static_cast<dd4hep::xml::Component>(xLayerColl);
    dd4hep::xml::Component xRods = xLayer.child("rods");
    dd4hep::xml::Component xRodEven = xRods.child("rodOdd");
    dd4hep::xml::Component xRodOdd = xRods.child("rodEven");
    dd4hep::xml::Component xModulesEven = xRodEven.child("modules");
    dd4hep::xml::Component xModulePropertiesOdd =
        xRodOdd.child("moduleProperties");
    dd4hep::xml::Component xModulesOdd = xRodOdd.child("modules");
    dd4hep::Tube layerShape(xLayer.rmin(), xLayer.rmax(), dimensions.zmax());
    Volume layerVolume("layer", layerShape, lcdd.material("Air"));
    // layerVolume.setVisAttributes(lcdd.invisible());
    PlacedVolume placedLayerVolume = topVolume.placeVolume(layerVolume);
    placedLayerVolume.addPhysVolID("layer", layerCounter);
    DetElement lay_det(topDetElement, "layer" + std::to_string(layerCounter),
                       layerCounter);
    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    // detElement owns extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive cylinder", "layer");
    layerExtension->addType("axes", "definitions", "XzY");
    lay_det.addExtension<Acts::ActsExtension>(layerExtension);
    lay_det.setPlacement(placedLayerVolume);
    dd4hep::xml::Component xModuleComponentsOdd =
        xModulePropertiesOdd.child("components");
    integratedModuleComponentThickness = 0;
    int moduleCounter = 0;
    Volume moduleVolume;

    // collect tracker material
    std::vector<std::pair<dd4hep::Material, double>> compMaterials;

    for (dd4hep::xml::Collection_t xModuleComponentOddColl(xModuleComponentsOdd,
                                                           _U(component));
         nullptr != xModuleComponentOddColl; ++xModuleComponentOddColl) {
      dd4hep::xml::Component xModuleComponentOdd =
          static_cast<dd4hep::xml::Component>(xModuleComponentOddColl);
      // collect module materials
      compMaterials.push_back(
          std::make_pair(lcdd.material(xModuleComponentOdd.materialStr()),
                         xModuleComponentOdd.thickness()));
    }
    for (dd4hep::xml::Collection_t xModuleComponentOddColl(xModuleComponentsOdd,
                                                           _U(component));
         nullptr != xModuleComponentOddColl; ++xModuleComponentOddColl) {
      dd4hep::xml::Component xModuleComponentOdd =
          static_cast<dd4hep::xml::Component>(xModuleComponentOddColl);
      auto moduleWidth = 0.5 * xModulePropertiesOdd.attr<double>("modWidth");
      auto moduleThickness = 0.5 * xModuleComponentOdd.thickness();
      auto moduleLength = 0.5 * xModulePropertiesOdd.attr<double>("modLength");
      Volume moduleVolume = Volume(
          "module", dd4hep::Box(moduleWidth, moduleThickness, moduleLength),
          lcdd.material(xModuleComponentOdd.materialStr()));

      // Create digitization module
      // with readout given by layer
      auto digiModule = det::utils::rectangleDigiModuleXZ(
          moduleWidth, moduleLength, moduleThickness, xLayer.X(), xLayer.Z());

      moduleVolume.setVisAttributes(lcdd.invisible());
      unsigned int nPhi = xRods.repeat();
      dd4hep::xml::Handle_t currentComp;
      for (unsigned int phiIndex = 0; phiIndex < nPhi; ++phiIndex) {
        double lX = 0;
        double lY = 0;
        double lZ = 0;
        if (0 == phiIndex % 2) {
          phi = 2 * M_PI * static_cast<double>(phiIndex) /
                static_cast<double>(nPhi);
          currentComp = xModulesEven;
        } else {
          currentComp = xModulesOdd;
        }
        for (dd4hep::xml::Collection_t xModuleColl(currentComp, _U(module));
             nullptr != xModuleColl; ++xModuleColl) {
          dd4hep::xml::Component xModule =
              static_cast<dd4hep::xml::Component>(xModuleColl);
          double currentPhi = atan2(xModule.Y(), xModule.X());
          double componentOffset =
              integratedModuleComponentThickness -
              0.5 * xModulePropertiesOdd.attr<double>("modThickness") +
              0.5 * xModuleComponentOdd.thickness();
          lX = xModule.X() + cos(currentPhi) * componentOffset;
          lY = xModule.Y() + sin(currentPhi) * componentOffset;
          lZ = xModule.Z();
          dd4hep::Translation3D moduleOffset(lX, lY, lZ);
          dd4hep::Transform3D lTrafo(
              dd4hep::RotationZ(atan2(lY, lX) + 0.5 * M_PI), moduleOffset);
          dd4hep::RotationZ lRotation(phi);
          PlacedVolume placedModuleVolume =
              layerVolume.placeVolume(moduleVolume, lRotation * lTrafo);
          if (xModuleComponentOdd.isSensitive()) {
            placedModuleVolume.addPhysVolID("module", moduleCounter);
            moduleVolume.setSensitiveDetector(sensDet);
            DetElement mod_det(lay_det,
                               "module" + std::to_string(moduleCounter),
                               moduleCounter);
            // add extension to hand over material
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
            mod_det.addExtension<Acts::ActsExtension>(moduleExtension);

            mod_det.setPlacement(placedModuleVolume);
            ++moduleCounter;
          }
        }
      }
      integratedModuleComponentThickness += xModuleComponentOdd.thickness();
    }
    ++layerCounter;
  }
  Volume motherVol = lcdd.pickMotherVolume(topDetElement);
  PlacedVolume placedGenericTrackerBarrel = motherVol.placeVolume(topVolume);
  placedGenericTrackerBarrel.addPhysVolID("system", topDetElement.id());
  topDetElement.setPlacement(placedGenericTrackerBarrel);
  return topDetElement;
}
}  // namespace det

DECLARE_DETELEMENT(TkLayoutBrlTracker, det::createTkLayoutTrackerBarrel)
