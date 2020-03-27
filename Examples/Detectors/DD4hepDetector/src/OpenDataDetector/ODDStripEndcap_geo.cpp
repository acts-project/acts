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

  Assembly diskAssembly("Disk");

  // DetElement tree
  DetElement diskElementTemplate("DiskElementTemplate", 0);

  // build the ring templates
  size_t ringNum = 0;
  for (xml_coll_t ring(x_det, _U(ring)); ring; ++ring, ++ringNum) {
    xml_comp_t x_ring = ring;

    string ringName = "Ring" + std::to_string(ringNum);
    Assembly ringAssembly(ringName);

    // DetElement tree
    DetElement ringElement(ringName, ringNum);

    if (x_ring.hasChild(_U(module))) {
      xml_comp_t x_module = x_ring.child(_U(module));
      auto module = assembleTrapezoidalModule(oddd, sens, x_module);

      // place the modules
      unsigned int nPhi = x_ring.nphi();
      double phiStep = 2 * M_PI / nPhi;
      double phi0 = x_ring.phi0();
      double r = x_ring.r();
      double zgap = x_ring.gap();
      bool reflect = x_ring.reflect();

      for (unsigned int modNum = 0; modNum < nPhi; ++modNum) {
        // The module name
        string moduleName = _toString((int)modNum, "module%d");

        bool odd = bool(modNum % 2);

        // Position parameters
        double phi = phi0 + modNum * phiStep;
        double x = r * cos(phi);
        double y = r * sin(phi);
        double z = odd ? -zgap : zgap;

        // Place Module Box Volumes, flip if necessary
        Position trans(x, y, z);
        double flip = odd ? M_PI : 0.;

        double angX = 0.5 * M_PI + flip;
        double angY = odd ? 0.5 * M_PI - phi : 0.5 * M_PI + phi;

        PlacedVolume placedModule = ringAssembly.placeVolume(
            module.first,
            Transform3D(
                RotationX(angX) * RotationY(angY),
                trans));  // RotationZ(phi + 1.5 * M_PI) * RotationY(flip)
        placedModule.addPhysVolID("module", modNum);
        // Clone the detector element
        auto moduleElement = module.second.clone(moduleName, modNum);
        moduleElement.setPlacement(placedModule);
        // Assign it as child to the stave template
        ringElement.add(moduleElement);
      }

      // Now add the ring detector Element to the disk
      diskElementTemplate.add(ringElement);

      size_t supportNum = 0;
      for (xml_coll_t sup(x_ring, _U(support)); sup; ++sup, ++supportNum) {
        xml_comp_t x_support = sup;
        // Create the volume of the support structure
        string supportName = _toString((int)supportNum, "RingSupport%d");
        Volume supportVolume(
            supportName,
            Tube(x_support.rmin(), x_support.rmax(), x_support.dz()),
            oddd.material(x_support.materialStr()));
        supportVolume.setVisAttributes(oddd, x_support.visStr());
        // Place the support structure
        PlacedVolume placedSupport = ringAssembly.placeVolume(
            supportVolume, Position(0., 0., x_support.z_offset()));
      }

      // Cooling rings
      buildCoolingRings(oddd, ringAssembly, x_ring);

      PlacedVolume placedRing = diskAssembly.placeVolume(
          ringAssembly, Position(0., 0., x_ring.z_offset()));

      placedRing.addPhysVolID("ring", ringNum);
      ringElement.setPlacement(placedRing);
    }
  }

  // Loop over the layers and place the disk, remember the z positions
  std::vector<double> endcapZ;
  size_t layNum = 0;
  for (xml_coll_t lay(xml, _U(layer)); lay; ++lay, ++layNum) {
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
    PlacedVolume placedLayer =
        endcapVolume.placeVolume(layerVolume, Position(0., 0., zeff));
    placedLayer.addPhysVolID("layer", layNum);
    endcapZ.push_back(zeff);

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
    }
    // Finish up the DetElement tree
    layerElement.setPlacement(placedLayer);
    endcapDetector.add(layerElement);
  }

  if (x_det.hasChild(_Unicode(services))) {
    // Grab the services - cables
    xml_comp_t x_services = x_det.child(_Unicode(services));
    for (xml_coll_t crout(x_services, _Unicode(cable_routing)); crout;
         ++crout) {
      xml_comp_t x_cable_routing = crout;
      buildEndcapRouting(oddd, endcapVolume, x_cable_routing, endcapZ);
    }
    // Grab for services - cooling
    for (xml_coll_t crout(x_services, _Unicode(cooling_routing)); crout;
         ++crout) {
      xml_comp_t x_cooling_routing = crout;
      buildEndcapRouting(oddd, endcapVolume, x_cooling_routing, endcapZ);
    }
  }

  // Place Volume
  Volume motherVolume = oddd.pickMotherVolume(endcapDetector);
  Position translation(0., 0., x_det_dim.z());
  PlacedVolume placedEndcap =
      motherVolume.placeVolume(endcapVolume, translation);
  // "system" is hard coded in the DD4Hep::VolumeManager
  placedEndcap.addPhysVolID("system", endcapDetector.id());
  endcapDetector.setPlacement(placedEndcap);

  // And return it
  return endcapDetector;
}

DECLARE_DETELEMENT(ODDStripEndcap, create_element)
