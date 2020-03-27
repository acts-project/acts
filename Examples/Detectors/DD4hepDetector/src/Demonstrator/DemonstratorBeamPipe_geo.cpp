// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();

  // Make DetElement
  DetElement beamPipeElement(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* beamPipeExtension = new Acts::ActsExtension();
  beamPipeExtension->addType("beampipe", "layer");
  beamPipeElement.addExtension<Acts::ActsExtension>(beamPipeExtension);

  // Make Volume
  xml_comp_t x_det_def = x_det.child(_U(description));
  Tube tube_shape(x_det_def.rmin(), x_det_def.rmax(), x_det_def.dz());
  Volume tube_vol(det_name, tube_shape, lcdd.material(x_det_def.materialStr()));
  tube_vol.setVisAttributes(lcdd, x_det.visStr());

  // Place it in the mother
  Volume mother_vol = lcdd.pickMotherVolume(beamPipeElement);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol);
  placedTube.addPhysVolID("BeamTube", beamPipeElement.id());
  beamPipeElement.setPlacement(placedTube);

  return beamPipeElement;
}

DECLARE_DETELEMENT(DemonstratorBeamPipe, create_element)
