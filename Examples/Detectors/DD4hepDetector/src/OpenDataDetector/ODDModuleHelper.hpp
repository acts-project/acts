// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

/// This method assemples a trapezoidal module for the
/// strip detectors
///
/// @param odd is the top level detector
/// @param sens is the top level sensitive detector container
/// @param x_module is the xml component describing the module
///
/// It excpects `module_component` xml childs
///
// @return a pair for a template module assembly and detector element
static std::pair<Assembly, DetElement> assembleTrapezoidalModule(
    Detector& oddd, SensitiveDetector& sens, const xml_comp_t& x_module) {
  // The Module envelope volume
  Assembly moduleAssembly("module");
  // Visualization
  moduleAssembly.setVisAttributes(oddd, x_module.visStr());

  // The module detector element
  DetElement moduleElement("ModuleElementTemplate", 0);

  // Place the components inside the module
  unsigned int compNum = 0;
  unsigned int sensorNum = 0;

  for (xml_coll_t comp(x_module, _U(module_component)); comp;
       ++comp, ++compNum) {
    xml_comp_t x_comp = comp;

    // create the component volume
    string compName =
        _toString((int)compNum, "component%d") + x_comp.materialStr();

    Trapezoid trapShape(x_comp.x1(), x_comp.x2(), 0.5 * x_comp.thickness(),
                        0.5 * x_comp.thickness(), x_comp.length());

    Volume componentVolume(compName, trapShape,
                           oddd.material(x_comp.materialStr()));
    componentVolume.setVisAttributes(oddd, x_comp.visStr());

    // overwrite if you have a subtraction
    if (x_comp.hasChild(_U(subtraction))) {
      xml_comp_t x_sub = x_comp.child(_U(subtraction));
      Tube tubeCutout(x_sub.rmin(), x_sub.rmax(), 1.1 * x_comp.length());

      // Create the subtraction
      componentVolume = Volume(
          compName,
          SubtractionSolid(
              trapShape, tubeCutout,
              Position(x_sub.x_offset(), x_sub.y_offset(), x_sub.z_offset())),
          oddd.material(x_comp.materialStr()));

      // place a fitting pipe if available
      if (x_comp.hasChild(_U(tube))) {
        xml_comp_t x_pipe = x_comp.child(_U(tube));
        Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(), x_comp.length());
        // Create the subtraction
        Volume pipeVolume("CoolingPipe", coolingPipe,
                          oddd.material(x_pipe.materialStr()));
        pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

        PlacedVolume pacedPipe = componentVolume.placeVolume(
            pipeVolume,
            Position(x_pipe.x_offset(), x_pipe.y_offset(), x_pipe.z_offset()));
      }
    }

    // Place the component
    double stereoAlpha = x_comp.alpha();
    PlacedVolume placedComponent = moduleAssembly.placeVolume(
        componentVolume,
        Transform3D(
            RotationY(stereoAlpha),
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset())));

    // Deal with the sensitive sensor
    if (x_comp.isSensitive()) {
      componentVolume.setSensitiveDetector(sens);
      placedComponent.addPhysVolID("sensor", sensorNum++);

      // Create the sensor element and place it
      string sensorName = _toString((int)sensorNum, "sensor%d");
      DetElement sensorElement(moduleElement, sensorName, sensorNum);
      sensorElement.setPlacement(placedComponent);

      // Add the sensor extension
      Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
      sensorExtension->addType("sensor", "detector");
      sensorExtension->addType("axes", "definitions", "XZY");
      // Set the extension
      sensorElement.addExtension<Acts::ActsExtension>(sensorExtension);
    }
  }

  // return the module assembly
  return std::pair<Assembly, DetElement>(moduleAssembly, moduleElement);
}

/// This method assemples a rectangular module for the
/// pixel and strip detectors
///
/// @param odd is the top level detector
/// @param sens is the top level sensitive detector container
/// @param x_module is the xml component describing the module
/// @param ylength[in,out] is the maximal length in y of all components
///        to be used for stave building
///
/// It excpects `module_component` xml childs
///
// @return a pair for a template module assembly and detector element
static std::pair<Assembly, DetElement> assembleRectangularModule(
    Detector& oddd, SensitiveDetector& sens, const xml_comp_t& x_module,
    double& ylength) {
  // The Module envelope volume
  Assembly moduleAssembly("module");
  // Visualization
  moduleAssembly.setVisAttributes(oddd, x_module.visStr());

  // The module detector element
  DetElement moduleElement("ModuleElementTemplate", 0);

  // Place the components inside the module
  unsigned int compNum = 0;
  unsigned int sensorNum = 0;

  for (xml_coll_t comp(x_module, _U(module_component)); comp;
       ++comp, ++compNum) {
    xml_comp_t x_comp = comp;

    // Component volume
    string componentName = _toString((int)compNum, "component%d");
    Box boxShape(0.5 * x_comp.dx(), 0.5 * x_comp.dy(), 0.5 * x_comp.dz());
    // Standard component volume without cutout
    Volume componentVolume(componentName, boxShape,
                           oddd.material(x_comp.materialStr()));

    // overwrite if you have a subtraction
    if (x_comp.hasChild(_U(subtraction))) {
      xml_comp_t x_sub = x_comp.child(_U(subtraction));
      Tube tubeCutout(x_sub.rmin(), x_sub.rmax(), x_comp.dy());

      // Create the subtraction
      componentVolume =
          Volume(componentName,
                 SubtractionSolid(
                     boxShape, tubeCutout,
                     Transform3D(RotationX(0.5 * M_PI),
                                 Position(x_sub.x_offset(), x_sub.y_offset(),
                                          x_sub.z_offset()))),
                 oddd.material(x_comp.materialStr()));

      // place a fitting pipe if available
      if (x_comp.hasChild(_U(tube))) {
        xml_comp_t x_pipe = x_comp.child(_U(tube));
        Tube coolingPipe(x_pipe.rmin(), x_pipe.rmax(), 0.5 * x_comp.dy());
        // Create the subtraction
        Volume pipeVolume("CoolingPipe", coolingPipe,
                          oddd.material(x_pipe.materialStr()));
        pipeVolume.setVisAttributes(oddd, x_pipe.visStr());

        PlacedVolume pacedPipe = componentVolume.placeVolume(
            pipeVolume,
            Transform3D(RotationX(0.5 * M_PI),
                        Position(x_pipe.x_offset(), x_pipe.y_offset(),
                                 x_pipe.z_offset())));
      }
    }
    componentVolume.setVisAttributes(oddd, x_comp.visStr());

    // Calculate the module dimension
    double cylength =
        2. * abs(std::copysign(0.5 * x_comp.dy(), x_comp.y_offset()) +
                 x_comp.y_offset());
    ylength = cylength > ylength ? cylength : ylength;

    // Visualization
    componentVolume.setVisAttributes(oddd, x_comp.visStr());
    // Place Module Box Volumes in layer
    double stereoAlpha = x_comp.alpha();
    PlacedVolume placedComponent = moduleAssembly.placeVolume(
        componentVolume,
        Transform3D(
            RotationZ(stereoAlpha),
            Position(x_comp.x_offset(), x_comp.y_offset(), x_comp.z_offset())));

    // Deal with the sensitive sensor
    if (x_comp.isSensitive()) {
      componentVolume.setSensitiveDetector(sens);
      placedComponent.addPhysVolID("sensor", sensorNum++);

      // Create the sensor element and place it
      string sensorName = _toString((int)sensorNum, "sensor%d");
      DetElement sensorElement(moduleElement, sensorName, sensorNum);
      sensorElement.setPlacement(placedComponent);

      // Add the sensor extension
      Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
      sensorExtension->addType("sensor", "detector");
      sensorExtension->addType("axes", "definitions", "XYZ");
      // Set the extension
      sensorElement.addExtension<Acts::ActsExtension>(sensorExtension);
    }
  }

  // return the module assembly
  return std::pair<Assembly, DetElement>(moduleAssembly, moduleElement);
}
