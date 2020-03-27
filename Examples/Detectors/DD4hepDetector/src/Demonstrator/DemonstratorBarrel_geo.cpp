// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>
#include <vector>

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepMaterial.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;

/**
 Constructor for a cylindrical barrel volume, possibly containing layers and the
 layers possibly containing modules.
 */

static Ref_t create_element(Detector& lcdd, xml_h xml,
                            dd4hep::SensitiveDetector sens) {
  xml_det_t x_det = xml;
  std::string barrelName = x_det.nameStr();
  // Make dd4hep::DetElement
  dd4hep::DetElement barrelDetector(barrelName, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension* barrelExtension = new Acts::ActsExtension();
  barrelExtension->addType("barrel", "detector");
  barrelDetector.addExtension<Acts::ActsExtension>(barrelExtension);

  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  dd4hep::Tube barrelShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  dd4hep::Volume barrelVolume(barrelName, barrelShape,
                              lcdd.air());  // air at the moment change later
  barrelVolume.setVisAttributes(lcdd, x_det.visStr());

  // Loop over the layers and build them
  for (xml_coll_t j(xml, _U(layer)); j; ++j) {
    // Get the layer xml configuration
    xml_comp_t x_layer = j;
    double rmin = x_layer.rmin();
    double rmax = x_layer.rmax();
    unsigned int layerNum = x_layer.id();
    // Create Volume for Layer
    std::string layerName = barrelName + _toString((int)layerNum, "layer%d");
    dd4hep::Volume layerVolume(layerName, Tube(rmin, rmax, x_layer.dz()),
                               lcdd.material(x_layer.materialStr()));
    dd4hep::DetElement layerElement(barrelDetector, layerName, layerNum);
    // Visualization
    layerVolume.setVisAttributes(lcdd, x_layer.visStr());

    unsigned int supportNum = 0;
    // Place the support cylinder
    if (x_layer.hasChild(_U(support))) {
      xml_comp_t x_support = x_layer.child(_U(support));

      // Create the volume of the support structure
      dd4hep::Volume supportVolume(
          "SupportStructure",
          Tube(x_support.rmin(), x_support.rmax(), x_support.dz()),
          lcdd.material(x_support.materialStr()));
      supportVolume.setVisAttributes(lcdd, x_support.visStr());
      // Place the support structure
      dd4hep::PlacedVolume placedSupport =
          layerVolume.placeVolume(supportVolume);
      placedSupport.addPhysVolID("support", supportNum++);
    }

    // Construct the volume
    if (x_layer.hasChild(_U(module))) {
      xml_comp_t x_module = x_layer.child(_U(module));
      // create the module volume and its corresponing component volumes first
      dd4hep::Assembly moduleAssembly("module");
      // Visualization
      moduleAssembly.setVisAttributes(lcdd, x_module.visStr());
      if (x_module.isSensitive()) {
        moduleAssembly.setSensitiveDetector(sens);
      }

      xml_comp_t x_mod_placement = x_module.child(_Unicode(placements));
      unsigned int nphi = x_mod_placement.nphi();
      double phi0 = x_mod_placement.phi0();
      double phiTilt = x_mod_placement.phi_tilt();
      double r = x_mod_placement.r();
      double deltaPhi = 2 * M_PI / nphi;

      // Place the components inside the module
      unsigned int compNum = 1;

      std::vector<PlacedVolume> sensComponents;

      for (xml_coll_t comp(x_module, _U(module_component)); comp;
           ++comp, ++compNum) {
        xml_comp_t x_comp = comp;
        // Component volume
        std::string componentName = _toString((int)compNum, "component%d");
        dd4hep::Volume componentVolume(
            componentName,
            Box(0.5 * x_comp.dx(), 0.5 * x_comp.dy(), 0.5 * x_comp.dz()),
            lcdd.material(x_comp.materialStr()));
        if (x_comp.isSensitive()) {
          componentVolume.setSensitiveDetector(sens);
        }

        // Visualization
        componentVolume.setVisAttributes(lcdd, x_comp.visStr());
        // Place Module Box Volumes in layer
        dd4hep::PlacedVolume placedComponent = moduleAssembly.placeVolume(
            componentVolume,
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset()));
        placedComponent.addPhysVolID("component", compNum);
        // Remember the sensitive components of this module
        if (x_comp.isSensitive()) {
          sensComponents.push_back(placedComponent);
        }
      }

      // Add cooling pipe
      if (x_module.hasChild(_U(tubs))) {
        xml_comp_t x_tubs = x_module.child(_U(tubs));
        dd4hep::Volume pipeVolume(
            "CoolingPipe", Tube(x_tubs.rmin(), x_tubs.rmax(), x_tubs.length()),
            lcdd.material(x_tubs.materialStr()));
        pipeVolume.setVisAttributes(lcdd, x_tubs.visStr());
        // Place the cooling pipe into the module
        dd4hep::PlacedVolume placedPipe = moduleAssembly.placeVolume(
            pipeVolume,
            Transform3D(RotationX(0.5 * M_PI) * RotationY(0.5 * M_PI),
                        Position(x_tubs.x_offset(), x_tubs.y_offset(),
                                 x_tubs.z_offset())));
        placedPipe.addPhysVolID("support", supportNum++);
      }

      // Add mount
      if (x_module.hasChild(_U(anchor))) {
        xml_comp_t x_trd = x_module.child(_U(anchor));
        // create the two shapes first
        dd4hep::Trapezoid mountShape(x_trd.x1(), x_trd.x2(), x_trd.length(),
                                     x_trd.length(), x_trd.dz());

        dd4hep::Volume mountVolume("ModuleMount", mountShape,
                                   lcdd.material(x_trd.materialStr()));

        // Place the mount onto the module
        dd4hep::PlacedVolume placedMount = moduleAssembly.placeVolume(
            mountVolume,
            Transform3D(RotationZ(0.5 * M_PI),
                        Position(x_trd.x_offset(), x_trd.y_offset(),
                                 x_trd.z_offset())));
        placedMount.addPhysVolID("support", supportNum++);
      }

      // Add cable
      if (x_module.hasChild(_U(box))) {
        xml_comp_t x_cab = x_module.child(_U(box));
        dd4hep::Volume cableVolume(
            "Cable", Box(0.5 * x_cab.dx(), 0.5 * x_cab.dy(), 0.5 * x_cab.dz()),
            lcdd.material(x_cab.materialStr()));
        // Visualization
        cableVolume.setVisAttributes(lcdd, x_cab.visStr());
        // Place Module Box Volumes in layer
        dd4hep::PlacedVolume placedCable = moduleAssembly.placeVolume(
            cableVolume,
            Transform3D(RotationX(x_cab.alpha()),
                        Position(x_cab.x_offset(), x_cab.y_offset(),
                                 x_cab.z_offset())));
        placedCable.addPhysVolID("support", supportNum++);
      }

      // Place the modules
      for (int iphi = 0; iphi < nphi; ++iphi) {
        double phi = phi0 + iphi * deltaPhi;
        std::string moduleName = layerName + _toString((int)iphi, "module%d");
        Position trans(r * cos(phi), r * sin(phi), 0.);
        // Create detector element
        dd4hep::DetElement moduleElement(layerElement, moduleName, iphi);
        // Place the sensitive inside here
        unsigned int ccomp = 1;
        for (auto& sensComp : sensComponents) {
          dd4hep::DetElement componentElement(moduleElement, "component",
                                              ccomp++);
          componentElement.setPlacement(sensComp);
          // Add the sensor extension
          Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
          sensorExtension->addType("sensor", "detector");
          // Set the extension
          componentElement.addExtension<Acts::ActsExtension>(sensorExtension);
        }

        // Place Module Box Volumes in layer
        dd4hep::PlacedVolume placedModule = layerVolume.placeVolume(
            moduleAssembly,
            Transform3D(RotationY(0.5 * M_PI) * RotationX(-phi - phiTilt),
                        trans));
        placedModule.addPhysVolID("module", iphi + 1);

        // Assign module dd4hep::DetElement to the placed module volume
        moduleElement.setPlacement(placedModule);
      }
    }

    // Configure the ACTS extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("barrel", "layer");
    // Add the proto layer material
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      xmlToProtoSurfaceMaterial(x_layer_material, *layerExtension,
                                "layer_material");
    }
    layerElement.addExtension<Acts::ActsExtension>(layerExtension);

    // Place layer volume
    dd4hep::PlacedVolume placedLayer = barrelVolume.placeVolume(layerVolume);
    placedLayer.addPhysVolID("layer", layerNum);
    // Assign layer dd4hep::DetElement to layer volume
    layerElement.setPlacement(placedLayer);
  }

  // Place Volume
  dd4hep::Volume motherVolume = lcdd.pickMotherVolume(barrelDetector);
  dd4hep::PlacedVolume placedBarrel = motherVolume.placeVolume(barrelVolume);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedBarrel.addPhysVolID("system", barrelDetector.id());
  barrelDetector.setPlacement(placedBarrel);

  return barrelDetector;
}

DECLARE_DETELEMENT(DemonstratorBarrel, create_element)
