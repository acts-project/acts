// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/detail/DD4hepAttributeHelpers.hpp"
#include "Acts/Plugins/DD4hep/detail/DD4hepBinningHelper.hpp"
#include "Acts/Plugins/DD4hep/detail/DD4hepConversionHelpers.hpp"

#include <DD4hep/DetFactoryHelper.h>
#include <DD4hep/Objects.h>
#include <DDRec/DetectorData.h>
#include <XML/Utilities.h>

#include "DD4hepTestsHelper.hpp"

using namespace std;
using namespace dd4hep;

/// Standard create_cylinder(...) create a simple cylinder layer
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_cylinder_layer(Detector &dd, xml_h xml,
                                   SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_env = x_det.child(_Unicode(envelope));

  // Make DetElement
  DetElement cylinderLayerElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, cylinderLayerElement);

  // Layer parameters
  auto &layerParams =
      DD4hepTestsHelper::ensureExtension<dd4hep::rec::VariantParameters>(
          cylinderLayerElement);

  // Check if the cylinder has a surface binning instruction
  if (x_det.hasChild(_Unicode(surface_binning))) {
    xml_comp_t sfBinning = x_det.child(_Unicode(surface_binning));
    Acts::detail::decodeBinning(layerParams, sfBinning, "surface_binning",
                                {"z", "phi"});
  }

  // Create the envelope
  DetElement envelopeElement(detName + "_envelope", x_det.id());
  Tube envelopeShape(detName + "_shape", x_det_env.rmin(), x_det_env.rmax(),
                     x_det_env.dz());
  Volume envelopeVolume(detName, envelopeShape,
                        dd.material(x_det_env.materialStr()));
  envelopeVolume.setVisAttributes(dd, x_det.visStr());

  // Count the sensors
  unsigned int sensorID = 1u;
  xml_comp_t x_det_modules = x_det.child(_Unicode(modules));
  for (xml_coll_t bmodule(x_det_modules, _U(box)); bmodule; ++bmodule) {
    xml_comp_t x_det_box = bmodule;

    // Due to convention this causes an axis flip on x
    Box boxShape(0.5 * x_det_box.dx(), 0.5 * x_det_box.dy(),
                 0.5 * x_det_box.dz());

    // Set an orientation
    DetElement boxElement(detName + "_module" + std::to_string(sensorID),
                          sensorID);

    Volume boxVolume(detName, boxShape, dd.material(x_det_box.materialStr()));
    boxVolume.setVisAttributes(dd, x_det.visStr());
    boxVolume.setSensitiveDetector(sens);

    PlacedVolume placedBox = envelopeVolume.placeVolume(
        boxVolume, DD4hepTestsHelper::createTransform(x_det_box));

    placedBox.addPhysVolID("sensor", sensorID++);
    boxElement.setPlacement(placedBox);
    // Add the module elements
    cylinderLayerElement.add(boxElement);
  }

  // Place the envelope into mother
  Volume motherVolume = dd.pickMotherVolume(cylinderLayerElement);
  PlacedVolume placedEnvelope = motherVolume.placeVolume(
      envelopeVolume, DD4hepTestsHelper::createTransform(x_det_env));
  placedEnvelope.addPhysVolID("system", cylinderLayerElement.id());
  cylinderLayerElement.setPlacement(placedEnvelope);

  return cylinderLayerElement;
}

DECLARE_DETELEMENT(CylinderLayer, create_cylinder_layer)

/// Standard create_cylinder(...) create a simple disc layer
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_disc_layer(Detector &dd, xml_h xml,
                               SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  xml_comp_t x_det_env = x_det.child(_Unicode(envelope));

  // Make DetElement
  DetElement discLayerElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, discLayerElement);

  auto &layerParams =
      DD4hepTestsHelper::ensureExtension<dd4hep::rec::VariantParameters>(
          discLayerElement);

  // Check if the disk has a surface binning instruction
  if (x_det.hasChild(_Unicode(surface_binning))) {
    xml_comp_t sfBinning = x_det.child(_Unicode(surface_binning));
    Acts::detail::decodeBinning(layerParams, sfBinning, "surface_binning",
                                {"r", "phi"});
  }

  // Create the envelope
  DetElement envelopeElement(detName + "_envelope", x_det.id());
  Tube envelopeShape(detName + "_shape", x_det_env.rmin(), x_det_env.rmax(),
                     x_det_env.dz());
  Volume envelopeVolume(detName, envelopeShape,
                        dd.material(x_det_env.materialStr()));
  envelopeVolume.setVisAttributes(dd, x_det.visStr());

  unsigned int sensorID = 1u;

  xml_comp_t x_det_modules = x_det.child(_Unicode(modules));

  for (xml_coll_t bmodule(x_det_modules, _U(trap)); bmodule; ++bmodule) {
    xml_comp_t x_det_trap = bmodule;

    // Due to convention this causes an axis flip on x
    Trapezoid trapShape(0.5 * x_det_trap.dz(), 0.5 * x_det_trap.dz(),
                        x_det_trap.x1(), x_det_trap.x2(),
                        0.5 * x_det_trap.dy());

    // Set an orientation
    DetElement trapElement(detName + "_module" + std::to_string(sensorID),
                           sensorID);
    auto &params =
        DD4hepTestsHelper::ensureExtension<dd4hep::rec::VariantParameters>(
            trapElement);
    params.set<std::string>("axis_definitions", "YZ");

    Volume trapVolume(detName, trapShape,
                      dd.material(x_det_trap.materialStr()));
    trapVolume.setVisAttributes(dd, x_det.visStr());
    trapVolume.setSensitiveDetector(sens);

    PlacedVolume placedTrap = envelopeVolume.placeVolume(
        trapVolume, DD4hepTestsHelper::createTransform(x_det_trap));

    placedTrap.addPhysVolID("sensor", sensorID++);
    trapElement.setPlacement(placedTrap);
    // Add the module elements
    discLayerElement.add(trapElement);
  }

  // Add the envlepe Element to the discLayerEnvelop

  // Place the envelope into mother
  Volume motherVolume = dd.pickMotherVolume(discLayerElement);
  PlacedVolume placedEnvelope = motherVolume.placeVolume(
      envelopeVolume, DD4hepTestsHelper::createTransform(x_det_env));
  placedEnvelope.addPhysVolID("system", discLayerElement.id());
  discLayerElement.setPlacement(placedEnvelope);

  return discLayerElement;
}

DECLARE_DETELEMENT(DiscLayer, create_disc_layer)
