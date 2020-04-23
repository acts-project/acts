// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/XMLElements.h"

namespace det {

/**
Factory for a shape from multiple cylinders. Meant for material approximations.
Expected xml structure:
<detector type="SimpleLayeredCylinder" ...>
  <dimensions rmin="..." rmax="..." dz="..." z_offset="..."> <!-- envelope -->
  <layer rmin="..." rmax="..." dz="..." z_offset="..." material="...">
  ...
  <layer rmin="..." rmax="..." dz="..." z_offset="..." material="...">
</detector>
@author: Joschka Lingemann
*/

static dd4hep::Ref_t createSimpleLayeredCylinder(
    dd4hep::Detector& lcdd, dd4hep::xml::Handle_t xmlElement,
    dd4hep::SensitiveDetector sensDet) {
  dd4hep::xml::DetElement xmlDet =
      static_cast<dd4hep::xml::DetElement>(xmlElement);
  std::string name = xmlDet.nameStr();
  dd4hep::DetElement detElement(name, xmlDet.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* detvolume = new Acts::ActsExtension();
  detvolume->addType("barrel", "detector");
  detElement.addExtension<Acts::ActsExtension>(detvolume);
  // Create volume
  dd4hep::Volume experimentalHall = lcdd.pickMotherVolume(detElement);
  xml_comp_t dimensions(xmlDet.dimensions());
  dd4hep::Tube envelope(dimensions.rmin(), dimensions.rmax(), dimensions.dz());
  dd4hep::Volume envVolume(name, envelope,
                           lcdd.material(dimensions.materialStr()));

  // Create layer cylinders with their respective material, etc
  size_t layerIdx = 0;
  // loop through layers
  for (xml_coll_t layerIt(xmlDet, _U(layer)); layerIt; ++layerIt) {
    xml_comp_t layerDet = layerIt;
    dd4hep::Tube layerShape(layerDet.rmin(), layerDet.rmax(), layerDet.dz());
    std::string layerName = dd4hep::xml::_toString(layerIdx, "layer%d");
    dd4hep::Volume layerVolume(layerName, layerShape, lcdd.air());
    // Create layer detector element
    dd4hep::DetElement lay_det(detElement, layerName, layerIdx);
    // @todo use material string again layer.attr<std::string>("material"))

    // envVolume.placeVolume(layerVolume,
    // dd4hep::Transform3D(dd4hep::RotationZ(0.),
    // transLayer));
    if (layerDet.hasAttr(_U(vis))) {
      layerVolume.setVisAttributes(lcdd, layerDet.visStr());
    }
    if (layerDet.hasAttr(_U(sensitive))) {
      layerVolume.setSensitiveDetector(sensDet);
    }
    // Set Acts Extension
    Acts::ActsExtension* detlayer = new Acts::ActsExtension();
    detlayer->addType("passive cylinder", "layer");
    lay_det.addExtension<Acts::ActsExtension>(detlayer);

    // place the layer into the mother volume with a possible translation
    dd4hep::Position transLayer(0., 0., layerDet.z_offset());
    dd4hep::PlacedVolume placedLayerVolume =
        envVolume.placeVolume(layerVolume, transLayer);
    // set volume ID
    placedLayerVolume.addPhysVolID("layer", layerIdx);
    lay_det.setPlacement(placedLayerVolume);
    layerIdx++;
  }

  dd4hep::Position trans(0., 0., dimensions.z_offset());
  dd4hep::PlacedVolume envPhys = experimentalHall.placeVolume(
      envVolume, dd4hep::Transform3D(dd4hep::RotationZ(0.), trans));
  envPhys.addPhysVolID("system", detElement.id());
  detElement.setPlacement(envPhys);
  detElement.setVisAttributes(lcdd, xmlDet.visStr(), envVolume);

  return detElement;
}
}  // namespace det

DECLARE_DETELEMENT(SimpleLayeredCylinder, det::createSimpleLayeredCylinder)
