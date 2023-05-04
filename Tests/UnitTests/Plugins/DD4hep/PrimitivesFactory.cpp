// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/detail/DD4hepAttributeHelpers.hpp"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Objects.h"
#include "XML/Utilities.h"

using namespace std;
using namespace dd4hep;

Transform3D createTransform(const xml_comp_t &x_det_comp) {
  // Build the transform - center def
  double cx = Acts::detail::attrValueOr<double>(x_det_comp, "cx", 0.);
  double cy = Acts::detail::attrValueOr<double>(x_det_comp, "cy", 0.);
  double cz = Acts::detail::attrValueOr<double>(x_det_comp, "cz", 0.);

  double xx = Acts::detail::attrValueOr<double>(x_det_comp, "xx", 1.);
  double xy = Acts::detail::attrValueOr<double>(x_det_comp, "xy", 0.);
  double xz = Acts::detail::attrValueOr<double>(x_det_comp, "xz", 0.);

  double yx = Acts::detail::attrValueOr<double>(x_det_comp, "yx", 0.);
  double yy = Acts::detail::attrValueOr<double>(x_det_comp, "yy", 1.);
  double yz = Acts::detail::attrValueOr<double>(x_det_comp, "yz", 0.);

  Position xAxis(xx, xy, xz);
  Position yAxis(yx, yy, yz);
  Position zAxis = xAxis.Cross(yAxis);
  double zx = zAxis.X();
  double zy = zAxis.Y();
  double zz = zAxis.Z();

  // Create the transform
  return Transform3D(xx, xy, xz, cx, yx, yy, yz, cy, zx, zy, zz, cz);
}

/// Standard create_cylinder(...) create a simple cylinder
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_cylinder(Detector &dd, xml_h xml,
                             SensitiveDetector /*sens*/) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_tubs = x_det.child(_U(tubs));

  // Make DetElement
  DetElement cylinderElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, cylinderElement);

  string shapeName = x_det_tubs.nameStr();
  double phiMin = Acts::detail::attrValueOr<double>(x_det_tubs, "phimin", 0.);
  double phiMax =
      Acts::detail::attrValueOr<double>(x_det_tubs, "phimax", 2 * M_PI);
  Tube tubeShape(shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(),
                 x_det_tubs.dz(), phiMin, phiMax);
  Volume tubeVolume(detName, tubeShape, dd.material(x_det_tubs.materialStr()));
  tubeVolume.setVisAttributes(dd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = dd.pickMotherVolume(cylinderElement);
  PlacedVolume placedTube = motherVolume.placeVolume(tubeVolume);
  placedTube.addPhysVolID(detName, cylinderElement.id());
  cylinderElement.setPlacement(placedTube);

  // And return the element for further parsing
  return cylinderElement;
}

DECLARE_DETELEMENT(Cylinder, create_cylinder)

/// Standard create_disc(...) create a simple disc
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_disc(Detector &dd, xml_h xml, SensitiveDetector /*sens*/) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_tubs = x_det.child(_U(tubs));

  // Make DetElement
  DetElement discElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, discElement);

  string shapeName = x_det_tubs.nameStr();
  double phiMin = Acts::detail::attrValueOr<double>(x_det_tubs, "phimin", 0.);
  double phiMax =
      Acts::detail::attrValueOr<double>(x_det_tubs, "phimax", 2 * M_PI);

  Tube discShape = Tube(shapeName, x_det_tubs.rmin(), x_det_tubs.rmax(),
                        x_det_tubs.dz(), phiMin, phiMax);

  Volume discVolume(detName, discShape, dd.material(x_det_tubs.materialStr()));
  discVolume.setVisAttributes(dd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = dd.pickMotherVolume(discElement);
  PlacedVolume placedTube = motherVolume.placeVolume(discVolume);
  placedTube.addPhysVolID(detName, discElement.id());
  discElement.setPlacement(placedTube);

  // And return the element for further parsing
  return discElement;
}

DECLARE_DETELEMENT(Disc, create_disc)

/// Standard create_rectangle(...) create a simple rectangle
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_rectangle(Detector &dd, xml_h xml,
                              SensitiveDetector /*sens*/) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_box = x_det.child(_U(box));

  // Make DetElement
  DetElement rectElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, rectElement);

  string shapeName = x_det_box.nameStr();
  Box rectShape(shapeName, 0.5 * x_det_box.dx(), 0.5 * x_det_box.dy(),
                0.5 * x_det_box.dz());

  Volume rectVolume(detName, rectShape, dd.material(x_det_box.materialStr()));
  rectVolume.setVisAttributes(dd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = dd.pickMotherVolume(rectElement);
  PlacedVolume placedRect =
      motherVolume.placeVolume(rectVolume, createTransform(x_det_box));
  placedRect.addPhysVolID(detName, rectElement.id());
  rectElement.setPlacement(placedRect);

  // And return the element for further parsing
  return rectElement;
}

DECLARE_DETELEMENT(Rectangle, create_rectangle)

/// Standard create_trapezoid (...) create a simple trapezoid
///
/// @param dd the detector to which this is addedded
/// @param xml the input xml element
/// @param sens is ignored
///
/// @return a reference counted DetElement
static Ref_t create_trapezoid(Detector &dd, xml_h xml,
                              SensitiveDetector /*sens*/) {
  xml_det_t x_det = xml;
  string detName = x_det.nameStr();

  // Make Volume
  xml_comp_t x_det_trap = x_det.child(_U(trap));

  // Make DetElement
  DetElement trapElement(detName, x_det.id());
  dd4hep::xml::setDetectorTypeFlag(xml, trapElement);

  string shapeName = x_det_trap.nameStr();

  // Due to convention this causes an axis flip on x
  Trapezoid trapShape(x_det_trap.x1(), x_det_trap.x2(), 0.5 * x_det_trap.dz(),
                      0.5 * x_det_trap.dz(), 0.5 * x_det_trap.dy());

  Volume trapVolume(detName, trapShape, dd.material(x_det_trap.materialStr()));
  trapVolume.setVisAttributes(dd, x_det.visStr());

  // Place it in the mother
  Volume motherVolume = dd.pickMotherVolume(trapElement);
  PlacedVolume placedTrap =
      motherVolume.placeVolume(trapVolume, createTransform(x_det_trap));
  placedTrap.addPhysVolID(detName, trapElement.id());
  trapElement.setPlacement(placedTrap);

  // And return the element for further parsing
  return trapElement;
}

DECLARE_DETELEMENT(Trapezoid, create_trapezoid)
