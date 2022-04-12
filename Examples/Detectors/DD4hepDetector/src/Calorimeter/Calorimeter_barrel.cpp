// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsDD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"

#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

static Ref_t create_element(Detector& oddd, xml_h xml, SensitiveDetector) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make DetElement
  DetElement CaloBarrel(detName, x_det.id());

  // add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension* barrelExtension = new Acts::ActsExtension();
  barrelExtension->addType("barrel", "detector");
  CaloBarrel.addExtension<Acts::ActsExtension>(barrelExtension);

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  Tube barrelShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume barrelVolume(detName, barrelShape, oddd.air());
  barrelVolume.setVisAttributes(oddd, x_det.visStr());

  // Create the tubes with their respective material, etc
  size_t VolIdx = 0;
  // loop through layers
  for (xml_coll_t tubIt(x_det, _U(tubs)); tubIt; ++tubIt) {
    xml_comp_t VolDet = tubIt;
    dd4hep::Tube VolShape(VolDet.rmin(), VolDet.rmax(), VolDet.dz());
    std::string VolName = dd4hep::xml::_toString(VolDet.nameStr());
    dd4hep::Volume Volume(VolName, VolShape,
                          oddd.material(VolDet.materialStr()));
    // Create volume detector element
    dd4hep::DetElement vol_det(CaloBarrel, VolName, VolIdx);

    if (VolDet.hasAttr(_U(vis))) {
      Volume.setVisAttributes(oddd, VolDet.visStr());
    }

    // Set Acts Extension
    Acts::ActsExtension* detVolume = new Acts::ActsExtension();
    detVolume->addType("passive cylinder", "volume");
    vol_det.addExtension<Acts::ActsExtension>(detVolume);

    // place the layer into the mother volume with a possible translation
    dd4hep::Position transvol(0., 0., 0.);
    dd4hep::PlacedVolume placedVolume =
        barrelVolume.placeVolume(Volume, transvol);
    // set volume ID
    placedVolume.addPhysVolID(VolName, VolIdx);
    vol_det.setPlacement(placedVolume);
    VolIdx++;
  }

  // Place it in the mother
  Volume motherVolume = oddd.pickMotherVolume(CaloBarrel);
  PlacedVolume placedTube = motherVolume.placeVolume(barrelVolume);
  placedTube.addPhysVolID(detName, CaloBarrel.id());
  CaloBarrel.setPlacement(placedTube);

  // And return the element for further parsing
  return CaloBarrel;
}

DECLARE_DETELEMENT(Calorimeter, create_element)
