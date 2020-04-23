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
 Constructor for a disc like endcap volume, possibly containing layers and the
 layers possibly containing modules. Both endcaps, the positive and negative can
 be build with this constructor. Atlas like style
 */

static Ref_t create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string det_name = x_det.nameStr();
  // Make DetElement
  DetElement cylinderVolume(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* detvolume = new Acts::ActsExtension();
  cylinderVolume.addExtension<Acts::ActsExtension>(detvolume);
  // make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  Tube tube_shape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume tube_vol(det_name, tube_shape,
                  lcdd.air());  // air at the moment change later
  tube_vol.setVisAttributes(lcdd, x_det_dim.visStr());
  // go trough possible layers
  size_t layer_num = 0;
  // if it is the positive or negative endcap
  int sign = 1;
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
    int module_num_num = 0;
    // go trough possible modules
    if (x_layer.hasChild(_U(module))) {
      for (xml_coll_t i(x_layer, _U(module)); i; i++) {
        xml_comp_t x_module = i;
        int repeat = x_module.repeat();
        double deltaphi = 2. * M_PI / repeat;
        double radius = x_module.radius();
        double slicedz = x_module.dz();

        size_t module_num = 0;

        // Create the module volume
        Volume mod_vol(
            "module",
            Trapezoid(x_module.x1(), x_module.x2(), x_module.thickness(),
                      x_module.thickness(), x_module.length()),
            lcdd.material(x_module.materialStr()));
        mod_vol.setVisAttributes(lcdd, x_module.visStr());

        // create the Acts::DigitizationModule (needed to do geometric
        // digitization) for all modules which have digitization module
        auto digiModule = FW::DD4hep::trapezoidalDigiModule(
            x_module.x1(), x_module.x2(), x_module.length(),
            x_module.thickness(), sens.readout().segmentation());

        // the sensitive placed components to be used later to create the
        // DetElements
        std::vector<PlacedVolume> sensComponents;
        // the possible digitization module
        std::shared_ptr<const Acts::DigitizationModule> digiComponent = nullptr;
        // go through possible components
        int comp_num = 0;
        for (xml_coll_t comp(x_module, _U(module_component)); comp; comp++) {
          xml_comp_t x_comp = comp;
          // create the component volume
          string comp_name =
              _toString((int)comp_num, "component%d") + x_comp.materialStr();
          Volume comp_vol(
              comp_name,
              Trapezoid(x_comp.x1(), x_comp.x2(), x_comp.thickness(),
                        x_comp.thickness(), x_comp.length()),
              lcdd.material(x_comp.materialStr()));
          comp_vol.setVisAttributes(lcdd, x_comp.visStr());

          // create the Acts::DigitizationModule (needed to do geometric
          // digitization) for all modules which have the sdigitization
          // compoenent
          digiComponent = FW::DD4hep::trapezoidalDigiModule(
              x_comp.x1(), x_comp.x2(), x_comp.length(), x_comp.thickness(),
              sens.readout().segmentation());

          // Set Sensitive Volumes sensitive
          if (x_comp.isSensitive())
            comp_vol.setSensitiveDetector(sens);

          // place component in module
          Position translation(0., x_comp.z(), 0.);
          PlacedVolume placed_comp = mod_vol.placeVolume(comp_vol, translation);
          if (x_comp.isSensitive())
            sensComponents.push_back(placed_comp);
          placed_comp.addPhysVolID("component", module_num);
          ++comp_num;
        }

        // Place the Modules
        for (int k = 0; k < repeat; k++) {
          string zname = _toString((int)k, "z%d");

          double phi = deltaphi / dd4hep::rad * k;
          string module_name =
              zname + _toString((int)(repeat * module_num_num + module_num),
                                "module%d");
          Position trans(radius * cos(phi), radius * sin(phi), slicedz);
          // Create the module DetElement
          DetElement mod_det(lay_det, module_name,
                             repeat * module_num_num + module_num);
          // Set Sensitive Volumes sensitive
          if (x_module.isSensitive()) {
            mod_vol.setSensitiveDetector(sens);
            // Create and attach the extension for DD4Hep/Acts conversion
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
            mod_det.addExtension<Acts::ActsExtension>(moduleExtension);
          }

          int comp_num = 0;
          for (auto& sensComp : sensComponents) {
            // Create DetElement
            DetElement comp_det(mod_det, "component", comp_num);
            // Create and attach the extension
            Acts::ActsExtension* moduleExtension = new Acts::ActsExtension();
            comp_det.addExtension<Acts::ActsExtension>(moduleExtension);
            comp_det.setPlacement(sensComp);
            comp_num++;
          }
          Rotation3D rotation1(1., 0., 0., 0., 1., 0., 0., 0., 1.);
          // Place Module Box Volumes in layer
          Transform3D transf1(RotationX(0.5 * M_PI) *
                                  RotationY(phi + 0.5 * M_PI) *
                                  RotationZ(0.1 * M_PI),
                              trans);
          PlacedVolume placedmodule =
              layer_vol.placeVolume(mod_vol, rotation1 * transf1);
          placedmodule.addPhysVolID("module",
                                    repeat * module_num_num + module_num);
          // assign module DetElement to the placed module volume
          mod_det.setPlacement(placedmodule);
          ++module_num;
        }
        ++module_num_num;
      }
    }

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive disk", "layer");
    lay_det.addExtension<Acts::ActsExtension>(layerExtension);

    // Placed Layer Volume
    Position layer_pos(0., 0., x_layer.z());
    PlacedVolume placedLayer = tube_vol.placeVolume(layer_vol, layer_pos);
    placedLayer.addPhysVolID("layer", layer_num);
    lay_det.setPlacement(placedLayer);
    ++layer_num;
  }
  // Place Volume
  // if it is the negative endcap the normal vector needs to point into the
  // outside
  Position endcap_translation(0., 0., x_det_dim.z());
  Rotation3D rotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  if (x_det_dim.z() < 0.) {
    rotation.SetComponents(1., 0., 0., 0., -1., 0., 0., 0., -1.);
  }
  Transform3D endcap_transform(rotation, endcap_translation);
  Volume mother_vol = lcdd.pickMotherVolume(cylinderVolume);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol, endcap_transform);
  placedTube.addPhysVolID("system", cylinderVolume.id());
  cylinderVolume.setPlacement(placedTube);

  return cylinderVolume;
}

DECLARE_DETELEMENT(ACTS_Endcap, create_element)
