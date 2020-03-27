// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

static Ref_t create_element(Detector& oddd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_tubs = x_det.child(_U(tubs));

  // Make DetElement
  DetElement cylinderElement(detName, x_det.id());

  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* pcExtension = new Acts::ActsExtension();
  bool isBeamPipe = x_det.hasChild(_U(beampipe));
  pcExtension->addType("passive cylinder", "layer");
  if (isBeamPipe) {
    pcExtension->addType("beampipe", "layer");
  }
  // Add the proto layer material
  for (xml_coll_t lmat(x_det_tubs, _Unicode(layer_material)); lmat; ++lmat) {
    xml_comp_t x_layer_material = lmat;
    xmlToProtoSurfaceMaterial(x_layer_material, *pcExtension, "layer_material");
  }
  cylinderElement.addExtension<Acts::ActsExtension>(pcExtension);

  string shapeName = x_det_tubs.nameStr();
  Tube tubeShape(shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(),
                 x_det_tubs.dz());
  Volume tubeVolume(detName, tubeShape,
                    oddd.material(x_det_tubs.materialStr()));
  tubeVolume.setVisAttributes(oddd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = oddd.pickMotherVolume(cylinderElement);
  PlacedVolume placedTube = motherVolume.placeVolume(tubeVolume);
  placedTube.addPhysVolID(detName, cylinderElement.id());
  cylinderElement.setPlacement(placedTube);

  // And return the element for further parsing
  return cylinderElement;
}

DECLARE_DETELEMENT(ODDCylinder, create_element)
