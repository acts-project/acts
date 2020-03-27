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
 layers possibly containing modules.
 */

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();
  // Make DetElement
  DetElement cylinderVolume(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
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
    // Create Volume for Layer
    string layer_name = det_name + _toString((int)layer_num, "layer%d");
    Volume layer_vol(layer_name, Tube(l_rmin, l_rmax, x_layer.z()),
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

      // Place the Modules in z
      size_t module_num = 0;
      // create the module volume and its corresponing component volumes first
      Volume mod_vol(
          "module",
          Box(x_module.length(), x_module.width(), x_module.thickness()),
          lcdd.material(x_module.materialStr()));
      // Visualization
      mod_vol.setVisAttributes(lcdd, x_module.visStr());
      // the sensitive placed components to be used later to create the
      // DetElements
      std::vector<PlacedVolume> sensComponents;
      int comp_num = 0;
      // go through module components
      for (xml_coll_t comp(x_module, _U(module_component)); comp; ++comp) {
        string component_name = _toString((int)comp_num, "component%d");
        xml_comp_t x_component = comp;
        Volume comp_vol(component_name,
                        Box(x_component.length(), x_component.width(),
                            x_component.thickness()),
                        lcdd.material(x_component.materialStr()));
        comp_vol.setVisAttributes(lcdd, x_component.visStr());

        // make sensitive components sensitive
        if (x_component.isSensitive())
          comp_vol.setSensitiveDetector(sens);

        // Place Component in Module
        Position trans(x_component.x(), 0., x_component.z());
        PlacedVolume placedcomponent = mod_vol.placeVolume(comp_vol, trans);
        placedcomponent.addPhysVolID("component", comp_num);
        if (x_component.isSensitive())
          sensComponents.push_back(placedcomponent);
        comp_num++;
      }
      // add possible trapezoidal shape with hole for cooling pipe
      if (x_module.hasChild(_U(subtraction))) {
        xml_comp_t x_sub = x_module.child(_U(subtraction));
        xml_comp_t x_trd = x_sub.child(_U(trd));
        xml_comp_t x_tubs = x_sub.child(_U(tubs));
        string component_name = _toString((int)comp_num, "component%d");
        // create the two shapes first
        Trapezoid trap_shape(x_trd.x1(), x_trd.x2(), x_trd.length(),
                             x_trd.length(), x_trd.thickness());
        Tube tubs_shape(x_tubs.rmin(), x_tubs.rmax(), x_tubs.dz());
        // create the subtraction
        Volume sub_vol("subtraction_components",
                       SubtractionSolid(trap_shape, tubs_shape,
                                        Transform3D(RotationX(0.5 * M_PI))),
                       lcdd.material(x_sub.materialStr()));
        sub_vol.setVisAttributes(lcdd, x_sub.visStr());
        // Place the volume in the module
        PlacedVolume placedSub = mod_vol.placeVolume(
            sub_vol, Transform3D(RotationZ(0.5 * M_PI) * RotationY(M_PI),
                                 Position(0., 0., x_sub.z())));
        placedSub.addPhysVolID("component", comp_num);
        comp_num++;
      }
      // add posibble cooling pipe
      if (x_module.hasChild(_U(tubs))) {
        xml_comp_t x_tubs = x_module.child(_U(tubs));
        string component_name = _toString((int)comp_num, "component%d");
        Volume pipe_vol("CoolingPipe",
                        Tube(x_tubs.rmin(), x_tubs.rmax(), x_tubs.dz()),
                        lcdd.material(x_tubs.materialStr()));
        pipe_vol.setVisAttributes(lcdd, x_tubs.visStr());
        // Place the cooling pipe into the module
        PlacedVolume placedPipe = mod_vol.placeVolume(
            pipe_vol, Transform3D(RotationX(0.5 * M_PI) * RotationY(0.5 * M_PI),
                                  Position(0., 0., x_tubs.z())));
        placedPipe.addPhysVolID("component", comp_num);
        comp_num++;
      }
      // Now place the module
      for (int k = -zrepeat; k <= zrepeat; k++) {
        double r = (l_rmax + l_rmin) * 0.5;
        if (k % 2 == 0)
          r += dr;
        // Place the modules in phi
        for (int i = 0; i < repeat; ++i) {
          double phi = deltaphi / dd4hep::rad * i;
          string module_name =
              layer_name + _toString((int)module_num, "module%d");
          Position trans(r * cos(phi), r * sin(phi), k * dz);
          // create detector element
          DetElement mod_det(lay_det, module_name, module_num);
          // Set Sensitive Volumes sensitive
          if (x_module.isSensitive()) {
            mod_vol.setSensitiveDetector(sens);
          }
          int comp_num = 0;
          for (auto& sensComp : sensComponents) {
            DetElement component_det(mod_det, "component", comp_num);
            component_det.setPlacement(sensComp);
            comp_num++;
          }
          // Place Module Box Volumes in layer
          PlacedVolume placedmodule = layer_vol.placeVolume(
              mod_vol,
              Transform3D(RotationX(-0.5 * M_PI) * RotationZ(-0.5 * M_PI) *
                              RotationX(phi - 0.6 * M_PI),
                          trans));
          placedmodule.addPhysVolID("module", module_num);
          // assign module DetElement to the placed module volume
          mod_det.setPlacement(placedmodule);
          ++module_num;
        }
      }
    }

    if (x_layer.hasChild(_U(tubs))) {
      xml_comp_t x_support = x_layer.child(_U(tubs));
      double srmin = l_rmin + x_support.offset();
      double srmax = srmin + x_support.thickness();
      // create the volume of the support structure
      Volume support_vol("SupportStructure", Tube(srmin, srmax, x_support.dz()),
                         lcdd.material(x_support.materialStr()));
      support_vol.setVisAttributes(lcdd, x_support.visStr());
      // place the support structure
      PlacedVolume placedSupport = layer_vol.placeVolume(support_vol);
      //   placedSupport.addPhysVolID("support", 1);
    }

    // set granularity of layer material mapping and where material should be
    // mapped hand over modules to ACTS

    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive cylinder", "layer");
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

DECLARE_DETELEMENT(ACTS_IBLSimpleBarrel, create_element)
