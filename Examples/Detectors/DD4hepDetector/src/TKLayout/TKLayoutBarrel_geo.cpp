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
 Constructor for a cylindrical barrel volume, possibly containing layers and the
 layers possibly containing modules. CMS like style
*/

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();
  // Make DetElement
  DetElement cylinderVolume(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* detvolume = new Acts::ActsExtension();
  detvolume->addType("barrel", "detector");
  cylinderVolume.addExtension<Acts::ActsExtension>(detvolume);
  // make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  Tube tube_shape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume tube_vol(det_name, tube_shape,
                  lcdd.air());  // air at the moment change later
  tube_vol.setVisAttributes(lcdd, x_det_dim.visStr());
  // go trough possible layers
  size_t layer_num = 0;
  for (xml_coll_t j(xml, _U(layer)); j; ++j) {
    xml_comp_t x_layer = j;
    double l_rmin = x_layer.inner_r();
    double l_rmax = x_layer.outer_r();
    double l_length = x_layer.z();
    // Create Volume and DetElement for Layer
    string layer_name = det_name + _toString((int)layer_num, "layer%d");
    Volume layer_vol(layer_name, Tube(l_rmin, l_rmax, l_length),
                     lcdd.material(x_layer.materialStr()));
    DetElement lay_det(cylinderVolume, layer_name, layer_num);
    // Visualization
    layer_vol.setVisAttributes(lcdd, x_layer.visStr());
    // go trough possible modules
    if (x_layer.hasChild(_U(module))) {
      xml_comp_t x_module = x_layer.child(_U(module));
      int repeat = x_module.repeat();
      double deltaphi = 2. * M_PI / repeat;
      // slices in z
      xml_comp_t x_slice = x_layer.child(_U(slice));
      int zrepeat = x_slice.repeat();
      double dz = x_slice.z();
      double dr = x_slice.dr();
      double offsetrz = x_slice.dz();
      size_t module_num = 0;
      // Creat the module volume
      Volume mod_vol(
          "module",
          Box(x_module.length(), x_module.width(), x_module.thickness()),
          lcdd.material(x_module.materialStr()));
      // Place the Modules in z
      for (int k = -zrepeat; k <= zrepeat; k++) {
        double r = (l_rmax + l_rmin) * 0.5;
        string zname = _toString((int)k, "z%d");
        if (k % 2 == 0)
          r -= offsetrz;
        // Place the modules in phi
        for (int i = 0; i < repeat; ++i) {
          double radius = r;
          if (i % 2 == 0)
            radius -= dr;
          // Visualization
          mod_vol.setVisAttributes(lcdd, x_module.visStr());
          double phi = deltaphi / dd4hep::rad * i;
          string module_name = zname + _toString((int)i, "module%d");
          Position trans(radius * cos(phi), radius * sin(phi), k * dz);
          // Create the module DetElement
          DetElement mod_det(lay_det, module_name, module_num);
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
              Transform3D(RotationX(-0.5 * M_PI) * RotationZ(-0.5 * M_PI) *
                              RotationX(phi - 0.5 * M_PI),
                          trans));
          placedmodule.addPhysVolID("module", module_num);
          // assign module DetElement to the placed module volume
          mod_det.setPlacement(placedmodule);
          ++module_num;
        }
      }
    }
    // set granularity of layer material mapping and where material should be
    // mapped
    // hand over modules to ACTS
    Acts::ActsExtension* detlayer = new Acts::ActsExtension();
    detlayer->addType("sensitive cylinder", "layer");
    lay_det.addExtension<Acts::ActsExtension>(detlayer);
    // Place layer volume
    PlacedVolume placedLayer = tube_vol.placeVolume(layer_vol);
    placedLayer.addPhysVolID("layer", layer_num);
    // Assign layer DetElement to layer volume
    lay_det.setPlacement(placedLayer);
    ++layer_num;
  }
  // Place Volume
  Volume mother_vol = lcdd.pickMotherVolume(cylinderVolume);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol);
  placedTube.addPhysVolID("system", cylinderVolume.id());
  cylinderVolume.setPlacement(placedTube);

  return cylinderVolume;
}

DECLARE_DETELEMENT(ACTS_TKLayoutBarrel, create_element)
