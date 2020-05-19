// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

/**
 Constructor for a disc like endcap volume, possibly containing layers and the
 layers possibly containing modules. Both endcaps, the positive and negative can
 be build with this constructor. CMS like style
 */

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();
  // Make DetElement
  DetElement cylinderVolume(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* detvolume = new Acts::ActsExtension();
  detvolume->addType("endcap", "detector");
  cylinderVolume.addExtension<Acts::ActsExtension>(detvolume);
  // make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  Tube tube_shape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume tube_vol(det_name, tube_shape,
                  lcdd.air());  // air at the moment change later
  tube_vol.setVisAttributes(lcdd, x_det_dim.visStr());
  // go trough possible layers
  int module_num_num = 0;
  size_t layer_num = 0;
  for (xml_coll_t j(xml, _U(layer)); j; ++j) {
    xml_comp_t x_layer = j;
    double l_rmin = x_layer.inner_r();
    double l_rmax = x_layer.outer_r();
    double l_length = x_layer.dz();
    // Create Volume and DetElement for Layer
    string layer_name = det_name + _toString((int)layer_num, "layer%d");
    Volume layer_vol(layer_name, Tube(l_rmin, l_rmax, l_length),
                     lcdd.material(x_layer.materialStr()));
    DetElement lay_det(cylinderVolume, layer_name, layer_num);
    // Visualization
    layer_vol.setVisAttributes(lcdd, x_layer.visStr());
    // go trough possible modules
    if (x_layer.hasChild(_U(module))) {
      for (xml_coll_t i(x_layer, _U(module)); i; i++) {
        xml_comp_t x_module = i;
        int repeat = x_module.repeat();
        double deltaphi = 2. * M_PI / repeat;
        double radius = x_module.radius();
        // Create the module volume
        Volume mod_vol(
            "module",
            Trapezoid(x_module.x1(), x_module.x2(), x_module.thickness(),
                      x_module.thickness(), x_module.length()),
            lcdd.material(x_module.materialStr()));
        size_t module_num = 0;
        // Place the Modules
        for (int k = 0; k < repeat; k++) {
          double slicedz = x_module.dz();
          if (k % 2 == 0)
            slicedz -= 10. * x_module.thickness();
          string zname = _toString((int)k, "z%d");
          // Visualization
          mod_vol.setVisAttributes(lcdd, x_module.visStr());
          double phi = deltaphi / dd4hep::rad * k;
          string module_name =
              zname + _toString((int)(repeat * module_num_num + module_num),
                                "module%d");
          Position trans(radius * cos(phi), radius * sin(phi), slicedz);
          // Create the module DetElement
          DetElement mod_det(lay_det, module_name,
                             repeat * module_num_num + module_num);
          // Create and attach the extension for DD4Hep/Acts conversion
          Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
          mod_det.addExtension<Acts::ActsExtension>(moduleExtension);
          // Set Sensitive Volmes sensitive
          if (x_module.isSensitive()) {
            mod_vol.setSensitiveDetector(sens);
          }
          // Place Module Box Volumes in layer
          PlacedVolume placedmodule = layer_vol.placeVolume(
              mod_vol,
              Transform3D(RotationX(0.5 * M_PI) * RotationY(phi + 0.5 * M_PI),
                          trans));
          placedmodule.addPhysVolID("module",
                                    repeat * module_num_num + module_num);
          // assign module DetElement to the placed module volume
          mod_det.setPlacement(placedmodule);
          ++module_num;
        }
        ++module_num_num;
      }
    }
    // set granularity of layer material mapping and where material should be
    // mapped
    // hand over modules to ACTS
    Acts::ActsExtension* detlayer = new Acts::ActsExtension();
    detlayer->addType("axes", "definitions", "XZy");
    detlayer->addType("sensitive disk", "layer");
    lay_det.addExtension<Acts::ActsExtension>(detlayer);
    double layerZpos = x_layer.z();
    // Placed Layer Volume
    Position layer_pos(0., 0., layerZpos);
    PlacedVolume placedLayer = tube_vol.placeVolume(layer_vol, layer_pos);
    placedLayer.addPhysVolID("layer", layer_num);
    lay_det.setPlacement(placedLayer);
    ++layer_num;
  }

  // if it is the negative endcap the normal vector needs to point into the
  // Place Volume
  Position endcap_translation(0., 0., x_det_dim.z());
  Rotation3D rotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  if (x_det_dim.z() < 0.) {
    rotation.SetComponents(1., 0., 0., 0., 1., 0., 0., 0., -1.);
  }
  Transform3D endcap_transform(rotation, endcap_translation);
  Volume mother_vol = lcdd.pickMotherVolume(cylinderVolume);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol, endcap_transform);
  placedTube.addPhysVolID("system", cylinderVolume.id());
  cylinderVolume.setPlacement(placedTube);

  return cylinderVolume;
}

DECLARE_DETELEMENT(ACTS_TKLayoutEndcap, create_element)
