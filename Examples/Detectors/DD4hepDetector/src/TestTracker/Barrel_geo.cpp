// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/DD4hepDetector/DD4hepDetectorHelper.hpp"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

/**
 Constructor for a cylindrical barrel volume, possibly containing layers and the
 layers possibly containing modules. In Atlas style
*/

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();
  // Make DetElement
  DetElement cylinderVolume(det_name, x_det.id());
  // Add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension* barrelExtension = new Acts::ActsExtension();
  barrelExtension->addType("barrel", "detector");
  cylinderVolume.addExtension<Acts::ActsExtension>(barrelExtension);
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
      // Create the module volume
      Volume mod_vol(
          "module",
          Box(x_module.width(), x_module.length(), x_module.thickness()),
          lcdd.material(x_module.materialStr()));

      // create the Acts::DigitizationModule (needed to do geometric
      // digitization) for all modules which have the same segmentation
      auto digiModule = FW::DD4hep::rectangleDigiModule(
          x_module.length(), x_module.width(), x_module.thickness(),
          sens.readout().segmentation());

      // Visualization
      mod_vol.setVisAttributes(lcdd, x_module.visStr());
      size_t module_num = 0;
      // Place the Modules in z
      for (int k = -zrepeat; k <= zrepeat; k++) {
        double r = (l_rmax + l_rmin) * 0.5;
        string zname = _toString((int)k, "z%d");
        if (k % 2 == 0)
          r += dr;
        // Place the modules in phi
        for (int i = 0; i < repeat; ++i) {
          double phi = deltaphi / dd4hep::rad * i;
          string module_name = zname + _toString((int)i, "module%d");
          Position trans(r * cos(phi), r * sin(phi), k * dz);
          // Create the module DetElement
          DetElement mod_det(lay_det, module_name, module_num);
          // Set Sensitive Volmes sensitive
          if (x_module.isSensitive()) {
            mod_vol.setSensitiveDetector(sens);
            // Create and attach the extension for DD4Hep/Acts conversion
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
            mod_det.addExtension<Acts::ActsExtension>(moduleExtension);
          }
          // Place Module Box Volumes in layer
          PlacedVolume placedmodule = layer_vol.placeVolume(
              mod_vol,
              Transform3D(RotationX(0.5 * M_PI) * RotationY(phi - 0.6 * M_PI),
                          trans));
          placedmodule.addPhysVolID("module", module_num);
          // assign module DetElement to the placed module volume
          mod_det.setPlacement(placedmodule);
          ++module_num;
        }
      }
    }

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("active cylinder", "layer");
    layerExtension->addType("axes", "definitions", "YxZ");
    lay_det.addExtension<Acts::ActsExtension>(layerExtension);
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

DECLARE_DETELEMENT(ACTS_Barrel, create_element)
