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
#include "ODDModuleHelper.hpp"
#include "ODDServiceHelper.hpp"

using namespace std;
using namespace dd4hep;

static Ref_t create_element(Detector& oddd, xml_h xml, SensitiveDetector sens) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make DetElement
  DetElement endcapDetector(detName, x_det.id());

  // Add Extension to DetElement for the RecoGeometry
  Acts::ActsExtension* endcapExtension = new Acts::ActsExtension();
  endcapExtension->addType("endcap", "detector");
  // Add the volume boundary material if configured
  for (xml_coll_t bmat(x_det, _Unicode(boundary_material)); bmat; ++bmat) {
    xml_comp_t x_boundary_material = bmat;
    xmlToProtoSurfaceMaterial(x_boundary_material, *endcapExtension,
                              "boundary_material");
  }
  endcapDetector.addExtension<Acts::ActsExtension>(endcapExtension);
  // Make Volume
  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  string endcapShapeName = x_det_dim.nameStr();

  Tube endcapShape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.dz());
  Volume endcapVolume(detName, endcapShape, oddd.air());
  endcapVolume.setVisAttributes(oddd, x_det.visStr());

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(endcapDetector);
  Position translation(0., 0., x_det_dim.z());
  PlacedVolume placedEndcap =
      motherVolume.placeVolume(endcapVolume, translation);

  // Create the module components
  xml_comp_t x_module = x_det.child(_U(module));
  double ylength = 0.;
  auto module = assembleRectangularModule(oddd, sens, x_module, ylength);

  Assembly diskAssembly("disk");

  // DetElement tree
  DetElement diskElementTemplate("DiskElementTemplate", 0);

  // Loop over the rings to create a template disk
  size_t ringNum = 0;
  for (xml_coll_t ring(xml, _U(ring)); ring; ++ring, ++ringNum) {
    // Get the ring
    xml_comp_t x_ring = ring;

    // The ring name
    string ringName = _toString((int)ringNum, "ring%d");
    Assembly ringAssembly(ringName);
    ringAssembly.setVisAttributes(oddd, x_ring.visStr());

    // DetElement tree
    DetElement ringElement(ringName, ringNum);

    double r = x_ring.r();
    double phi0 = x_ring.phi0();
    unsigned int nModules = x_ring.nphi();
    double zgap = x_ring.gap();
    double phiStep = 2. * M_PI / nModules;

    // Loop over modules
    for (unsigned int modNum = 0; modNum < nModules; ++modNum) {
      // The module name
      string moduleName = _toString((int)modNum, "module%d");
      bool odd = bool(modNum % 2);
      // Position it
      double phi = phi0 + modNum * phiStep;
      double z = odd ? -zgap : zgap;
      Position trans(r * cos(phi), r * sin(phi), z);
      // Place Module Box Volumes, flip if necessary
      double flip = x_det_dim.z() < 0. ? M_PI : 0.;
      if (ringNum != 0) {
        flip += M_PI;
      }
      PlacedVolume placedModule = ringAssembly.placeVolume(
          module.first,
          Transform3D(RotationZ(phi + 1.5 * M_PI) * RotationY(flip), trans));
      placedModule.addPhysVolID("module", modNum);
      // Clone the detector element
      auto moduleElement = module.second.clone(moduleName, modNum);
      moduleElement.setPlacement(placedModule);
      // Assign it as child to the stave template
      ringElement.add(moduleElement);
    }

    // Place Ring assembly into disk
    PlacedVolume placedRing = diskAssembly.placeVolume(
        ringAssembly, Position(0., 0., x_ring.z_offset()));
    placedRing.addPhysVolID("ring", ringNum);
    ringElement.setPlacement(placedRing);
    // Add it to the Disk element template
    diskElementTemplate.add(ringElement);
  }

  xml_comp_t x_support = x_det.child(_Unicode(ring_support));
  // The support shape
  Tube supportShape(x_support.rmin(), x_support.rmax(), x_support.dz());
  Volume supportVolume("DiskSupport", supportShape,
                       oddd.material(x_support.materialStr()));
  supportVolume.setVisAttributes(oddd, x_support.visStr());
  diskAssembly.placeVolume(supportVolume);

  // Cooling rings
  buildCoolingRings(oddd, diskAssembly, x_det);

  // Loop over the layers and place the disk
  size_t layNum = 0;
  // Remember the layers for the service routing
  std::vector<double> endcapZ;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layNum) {
    // Get the layer
    xml_comp_t x_layer = lay;

    // The Layer envelope volume
    string layerName = detName + std::to_string(layNum);
    Volume layerVolume(layerName,
                       Tube(x_layer.rmin(), x_layer.rmax(), x_layer.dz()),
                       oddd.air());

    layerVolume.setVisAttributes(oddd, x_layer.visStr());

    string diskElName = _toString((int)layNum, "disk%d");

    // The DetElement tree
    DetElement layerElement(layerName, layNum);
    auto diskElement = diskElementTemplate.clone(diskElName, layNum);

    // Place the disk into the layer
    PlacedVolume placedDisk = layerVolume.placeVolume(diskAssembly);
    diskElement.setPlacement(placedDisk);
    layerElement.add(diskElement);

    // Place Ring assembly into disk
    double zeff = x_layer.z_offset() - x_det_dim.z();
    endcapZ.push_back(zeff);

    PlacedVolume placedLayer =
        endcapVolume.placeVolume(layerVolume, Position(0., 0., zeff));
    placedLayer.addPhysVolID("layer", layNum);

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive disk", "layer");
    layerElement.addExtension<Acts::ActsExtension>(layerExtension);
    // Add the proto layer material
    for (xml_coll_t lmat(x_layer, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      xmlToProtoSurfaceMaterial(x_layer_material, *layerExtension,
                                "layer_material");
    }  // Finish up the DetElement tree
    layerElement.setPlacement(placedLayer);
    endcapDetector.add(layerElement);
  }

  // Close up the detector
  if (x_det.hasChild(_U(disk))) {
    // Endplate disk
    xml_comp_t x_endplate = x_det.child(_U(disk));

    // The Shape and Volume
    Tube endplateShape(x_endplate.rmin(), x_endplate.rmax(), x_endplate.dz());
    Volume endplateVolume("Endplate", endplateShape,
                          oddd.material(x_endplate.materialStr()));
    endplateVolume.setVisAttributes(oddd, x_endplate.visStr());

    double zeff = x_endplate.z_offset() - x_det_dim.z();
    endcapZ.push_back(zeff);
    PlacedVolume placedEndplate =
        endcapVolume.placeVolume(endplateVolume, Position(0., 0., zeff));

    DetElement endplateElement("Endplate", 0);

    // Place the layer with appropriate Acts::Extension
    // Configure the ACTS extension
    Acts::ActsExtension* endplateExtension = new Acts::ActsExtension();
    endplateExtension->addType("passive disk", "layer");
    endplateElement.addExtension<Acts::ActsExtension>(endplateExtension);
    // Add the proto layer material
    for (xml_coll_t lmat(x_endplate, _Unicode(layer_material)); lmat; ++lmat) {
      xml_comp_t x_layer_material = lmat;
      xmlToProtoSurfaceMaterial(x_layer_material, *endplateExtension,
                                "layer_material");
    }
    // Finish up the DetElement tree
    endplateElement.setPlacement(placedEndplate);
    endcapDetector.add(endplateElement);
  }

  if (x_det.hasChild(_Unicode(services))) {
    // Grab the services
    xml_comp_t x_services = x_det.child(_Unicode(services));
    if (x_services.hasChild(_Unicode(cable_routing))) {
      xml_comp_t x_cable_routing = x_services.child(_Unicode(cable_routing));
      buildEndcapRouting(oddd, endcapVolume, x_cable_routing, endcapZ);
    }
    if (x_services.hasChild(_Unicode(cooling_routing))) {
      xml_comp_t x_cooling_routing =
          x_services.child(_Unicode(cooling_routing));
      buildEndcapRouting(oddd, endcapVolume, x_cooling_routing, endcapZ);
    }
  }

  // Place the additional support cylinders per detector
  std::vector<double> layerR;
  buildSupportCylinder(oddd, endcapVolume, x_det, layerR);

  // "system" is hard coded in the DD4Hep::VolumeManager
  placedEndcap.addPhysVolID("system", endcapDetector.id());
  endcapDetector.setPlacement(placedEndcap);

  // And return it
  return endcapDetector;
}

DECLARE_DETELEMENT(ODDPixelEndcap, create_element)
