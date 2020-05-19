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
using dd4hep::xml::Component;
using dd4hep::xml::Dimension;

namespace det {
static dd4hep::Ref_t createTkLayoutTrackerEndcap(
    dd4hep::Detector& lcdd, dd4hep::xml::Handle_t xmlElement,
    dd4hep::SensitiveDetector sensDet) {
  // shorthands
  dd4hep::xml::DetElement xmlDet =
      static_cast<dd4hep::xml::DetElement>(xmlElement);
  Dimension dimensions(xmlDet.dimensions());
  double l_overlapMargin = 0.01;

  // get sensitive detector type from xml
  dd4hep::xml::Dimension sdTyp =
      xmlElement.child(_Unicode(sensitive));  // retrieve the type
  sensDet.setType(sdTyp.typeStr());           // set for the whole detector

  // definition of top volume
  std::string detName = xmlDet.nameStr();
  DetElement worldDetElement(detName, xmlDet.id());
  DetElement posEcapDetElement(worldDetElement, "posEndcap", 0);

  Acts::ActsExtension* ecapDetExt = new Acts::ActsExtension();
  ecapDetExt->addType("endcap", "detector");
  posEcapDetElement.addExtension<Acts::ActsExtension>(ecapDetExt);

  dd4hep::Assembly envelopeVolume("endcapEnvelope");
  envelopeVolume.setVisAttributes(lcdd.invisible());

  Component xDiscs = xmlElement.child(_Unicode(discs));

  double envelopeThickness = 0.5 * (dimensions.zmax() - dimensions.zmin());

  l_overlapMargin *= 0.9;

  unsigned int discCounter = 0;
  unsigned int compCounter = 0;
  double currentZ;
  std::vector<Volume> discVolumeVec;
  std::vector<DetElement> discDetElementVec;
  /// iterate over discs
  for (dd4hep::xml::Collection_t xDiscColl(xDiscs, _Unicode(discZPls));
       nullptr != xDiscColl; ++xDiscColl) {
    Component xDisc = static_cast<Component>(xDiscColl);
    Component xCurrentRings = xDisc.child(_Unicode(rings));
    // create disc volume
    double discThickness = 0.5 * (xDisc.zmax() - xDisc.zmin());
    currentZ = xDisc.z() - dimensions.zmin() - envelopeThickness;
    if (xCurrentRings.hasChild(
            _Unicode(ring))) {  // we have information to construct a new volume
      dd4hep::Tube discShape(xDisc.rmin() - l_overlapMargin,
                             xDisc.rmax() + l_overlapMargin,
                             discThickness + l_overlapMargin);

      discVolumeVec.emplace_back("disc", discShape, lcdd.air());
      discDetElementVec.emplace_back(
          posEcapDetElement, "disc" + std::to_string(discCounter), discCounter);
      // the local coordinate systems of modules in dd4hep and acts differ
      // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
      Acts::ActsExtension* detlayer = new Acts::ActsExtension();
      detlayer->addType("sensitive disk", "layer");
      detlayer->addType("axes", "definitions", "XZY");
      discDetElementVec.back().addExtension<Acts::ActsExtension>(detlayer);
      // iterate over rings
      for (dd4hep::xml::Collection_t xRingColl(xCurrentRings, _U(ring));
           (nullptr != xRingColl); ++xRingColl) {
        Component xRing = static_cast<Component>(xRingColl);
        Component xRingModules = xRing.child(_Unicode(modules));
        Component xModuleOdd = xRingModules.child(_Unicode(moduleOdd));
        Component xModuleEven = xRingModules.child(_Unicode(moduleEven));
        Component xModuleProperties = xRing.child(_Unicode(moduleProperties));
        Component xModulePropertiesComp =
            xModuleProperties.child(_Unicode(components));
        Component xSensorProperties = xRing.child(_Unicode(sensorProperties));

        // the component materials
        std::vector<std::pair<dd4hep::Material, double>> compMaterials;
        // place components in module
        for (dd4hep::xml::Collection_t xCompColl(xModulePropertiesComp,
                                                 _U(component));
             nullptr != xCompColl; ++xCompColl) {
          dd4hep::xml::Component xComp = static_cast<Component>(xCompColl);
          // collect module materials
          compMaterials.push_back(std::make_pair(
              lcdd.material(xComp.materialStr()), xComp.thickness()));
        }
        double integratedCompThickness = 0.;
        for (dd4hep::xml::Collection_t xCompColl(xModulePropertiesComp,
                                                 _U(component));
             nullptr != xCompColl; ++xCompColl) {
          Component xComp = static_cast<Component>(xCompColl);
          double compMinWidth =
              0.5 * xModuleProperties.attr<double>("modWidthMin");
          double compMaxWidth =
              0.5 * xModuleProperties.attr<double>("modWidthMax");
          double compThickness = 0.5 * xComp.thickness();
          double compLength =
              0.5 * xSensorProperties.attr<double>("sensorLength");
          Volume componentVolume(
              "component",
              dd4hep::Trapezoid(compMinWidth, compMaxWidth, compThickness,
                                compThickness, compLength),
              lcdd.material(xComp.materialStr()));

          // Create digitization module
          auto digiModule = det::utils::trapezoidalDigiModuleXZ(
              compMinWidth, compMaxWidth, compLength, compThickness, xRing.X(),
              xRing.Z());

          componentVolume.setVisAttributes(lcdd.invisible());
          unsigned int nPhi = xRing.attr<int>("nModules");
          double phi = 0;
          for (unsigned int phiIndex = 0; phiIndex < nPhi; ++phiIndex) {
            double lX = 0;
            double lY = 0;
            double lZ = 0;
            double phiTilt = 0;
            double thetaTilt = 0;
            if (0 == phiIndex % 2) {
              // the rotation for the odd module is already taken care
              // of by the position in tklayout xml
              phi = 2 * dd4hep::pi * static_cast<double>(phiIndex) /
                    static_cast<double>(nPhi);
              lX = xModuleEven.X();
              lY = xModuleEven.Y();
              lZ = xModuleEven.Z() - xDisc.zmin() - discThickness;
              phiTilt = xModuleEven.attr<double>("phiTilt");
              thetaTilt = xModuleEven.attr<double>("thetaTilt");
            } else {
              lX = xModuleOdd.X();
              lY = xModuleOdd.Y();
              lZ = xModuleOdd.Z() - xDisc.zmin() - discThickness;
              phiTilt = xModuleOdd.attr<double>("phiTilt");
              thetaTilt = xModuleOdd.attr<double>("thetaTilt");
            }
            // position module in the x-y plane, smaller end inward
            // and incorporate phi tilt if any
            dd4hep::RotationY lRotation1(M_PI * 0.5);
            dd4hep::RotationX lRotation2(M_PI * 0.5 + phiTilt);
            // align radially
            double componentOffset =
                integratedCompThickness -
                0.5 * xModuleProperties.attr<double>("modThickness") +
                0.5 * xComp.thickness();
            dd4hep::RotationZ lRotation3(atan2(lY, lX));
            // theta tilt, if any -- note the different convention between
            // tklayout and here, thus the subtraction of pi / 2
            dd4hep::RotationY lRotation4(thetaTilt - M_PI * 0.5);
            dd4hep::RotationZ lRotation_PhiPos(phi);
            // position in  disk
            dd4hep::Translation3D lTranslation(lX, lY, lZ + componentOffset);
            dd4hep::Transform3D myTrafo(
                lRotation4 * lRotation3 * lRotation2 * lRotation1,
                lTranslation);
            PlacedVolume placedComponentVolume =
                discVolumeVec.back().placeVolume(componentVolume,
                                                 lRotation_PhiPos * myTrafo);
            if (xComp.isSensitive()) {
              placedComponentVolume.addPhysVolID("component", compCounter);
              componentVolume.setSensitiveDetector(sensDet);
              DetElement moduleDetElement(discDetElementVec.back(),
                                          "comp" + std::to_string(compCounter),
                                          compCounter);

              // add extension to hand over material
              Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
              moduleDetElement.addExtension<Acts::ActsExtension>(
                  moduleExtension);

              moduleDetElement.setPlacement(placedComponentVolume);
              ++compCounter;
            }
          }
          integratedCompThickness += xComp.thickness();
        }
      }
    } else {
      discDetElementVec.emplace_back(
          discDetElementVec.back().clone("disc" + std::to_string(discCounter)));
      posEcapDetElement.add(discDetElementVec.back());
    }
    PlacedVolume placedDiscVolume = envelopeVolume.placeVolume(
        discVolumeVec.back(), dd4hep::Position(0, 0, currentZ));
    placedDiscVolume.addPhysVolID("disc", discCounter);
    ++discCounter;

    discDetElementVec.back().setPlacement(placedDiscVolume);
  }
  dd4hep::Assembly bothEndcaps("bothEndcapsEnvelope");

  dd4hep::Translation3D envelopeTranslation(
      0, 0, dimensions.zmin() + envelopeThickness);

  dd4hep::RotationX envelopeNegRotation(dd4hep::pi);
  dd4hep::RotationX envelopePosRotation(0.);
  PlacedVolume placedEnvelopeVolume = bothEndcaps.placeVolume(
      envelopeVolume, envelopePosRotation * envelopeTranslation);
  PlacedVolume placedNegEnvelopeVolume = bothEndcaps.placeVolume(
      envelopeVolume, envelopeNegRotation * envelopeTranslation);
  placedEnvelopeVolume.addPhysVolID("posneg", 0);
  placedNegEnvelopeVolume.addPhysVolID("posneg", 1);
  auto negEcapDetElement = posEcapDetElement.clone("negEndcap");

  posEcapDetElement.setPlacement(placedEnvelopeVolume);
  negEcapDetElement.setPlacement(placedNegEnvelopeVolume);
  worldDetElement.add(negEcapDetElement);
  // top of the hierarchy
  PlacedVolume mplv =
      lcdd.pickMotherVolume(worldDetElement).placeVolume(bothEndcaps);
  worldDetElement.setPlacement(mplv);
  mplv.addPhysVolID("system", xmlDet.id());
  return worldDetElement;
}
}  // namespace det

DECLARE_DETELEMENT(TkLayoutEcapTracker, det::createTkLayoutTrackerEndcap)
