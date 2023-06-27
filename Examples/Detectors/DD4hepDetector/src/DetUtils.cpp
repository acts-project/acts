// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "DetUtils.h"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cstddef>

#include <DD4hep/CartesianGridXZ.h>
#include <DD4hep/Segmentations.h>
#include <XML/XMLTags.h>

namespace det {
namespace utils {

std::shared_ptr<const Acts::DigitizationModule> rectangleDigiModuleXZ(
    double halflengthX, double halflengthZ, double thickness,
    const dd4hep::Segmentation& segmentation) {
  // convert to ACTS units
  double scalor = Acts::UnitConstants::cm;
  halflengthX *= scalor;
  halflengthZ *= scalor;
  thickness *= scalor;
  auto bounds =
      std::make_shared<const Acts::RectangleBounds>(halflengthX, halflengthZ);
  dd4hep::CartesianGridXZ cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    double gridSizeX = cartesianGrid.gridSizeX() * scalor;
    double gridSizeZ = cartesianGrid.gridSizeZ() * scalor;
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
                       ? static_cast<size_t>((2 * halflengthX) / gridSizeX)
                       : 0;
    size_t bins1 = (cartesianGrid.gridSizeZ() != 0)
                       ? static_cast<size_t>((2 * halflengthZ) / gridSizeZ)
                       : 0;

    std::shared_ptr<const Acts::CartesianSegmentation> actsSegmentation =
        std::make_shared<const Acts::CartesianSegmentation>(bounds, bins0,
                                                            bins1);
    // finally create the digitization module
    // @todo set lorentz angle
    return (std::make_shared<const Acts::DigitizationModule>(actsSegmentation,
                                                             thickness, 1, 0));
  }
  return nullptr;
}

std::shared_ptr<const Acts::DigitizationModule> rectangleDigiModuleXZ(
    double halflengthX, double halflengthZ, double thickness, double gridSizeX,
    double gridSizeZ) {
  // convert to ACTS units
  double scalor = Acts::UnitConstants::cm;
  halflengthX *= scalor;
  halflengthZ *= scalor;
  thickness *= scalor;
  auto bounds =
      std::make_shared<const Acts::RectangleBounds>(halflengthX, halflengthZ);
  // the Acts segmentation of the DigitizationModule
  size_t bins0 =
      (gridSizeX != 0)
          ? static_cast<size_t>((2 * halflengthX) / (gridSizeX * scalor))
          : 0;
  size_t bins1 =
      (gridSizeZ != 0)
          ? static_cast<size_t>((2 * halflengthZ) / (gridSizeZ * scalor))
          : 0;

  std::shared_ptr<const Acts::CartesianSegmentation> actsSegmentation =
      std::make_shared<const Acts::CartesianSegmentation>(bounds, bins0, bins1);

  // finally create the digitization module
  // @todo set lorentz angle
  return (std::make_shared<const Acts::DigitizationModule>(actsSegmentation,
                                                           thickness, 1, 0));
}

std::shared_ptr<const Acts::DigitizationModule> trapezoidalDigiModuleXZ(
    double minHalflengthX, double maxHalflengthX, double halflengthZ,
    double thickness, const dd4hep::Segmentation& segmentation) {
  // convert to ACTS units
  double scalor = Acts::UnitConstants::cm;
  minHalflengthX *= scalor;
  maxHalflengthX *= scalor;
  halflengthZ *= scalor;
  thickness *= scalor;

  auto bounds = std::make_shared<const Acts::TrapezoidBounds>(
      minHalflengthX, maxHalflengthX, halflengthZ);

  dd4hep::CartesianGridXZ cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    double gridSizeX = cartesianGrid.gridSizeX() * scalor;
    double gridSizeZ = cartesianGrid.gridSizeZ() * scalor;
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
                       ? static_cast<size_t>((2 * maxHalflengthX) / gridSizeX)
                       : 0;
    size_t bins1 = (cartesianGrid.gridSizeZ() != 0)
                       ? static_cast<size_t>((2 * halflengthZ) / gridSizeZ)
                       : 0;

    std::shared_ptr<const Acts::CartesianSegmentation> actsSegmentation =
        std::make_shared<const Acts::CartesianSegmentation>(bounds, bins0,
                                                            bins1);
    // finally create the digitization module
    // @todo set lorentz angle
    return (std::make_shared<const Acts::DigitizationModule>(actsSegmentation,
                                                             thickness, 1, 0));
  }
  return nullptr;
}

std::shared_ptr<const Acts::DigitizationModule> trapezoidalDigiModuleXZ(
    double minHalflengthX, double maxHalflengthX, double halflengthZ,
    double thickness, double gridSizeX, double gridSizeZ) {
  // convert to ACTS units
  double scalor = Acts::UnitConstants::cm;
  minHalflengthX *= scalor;
  maxHalflengthX *= scalor;
  halflengthZ *= scalor;
  thickness *= scalor;

  auto bounds = std::make_shared<const Acts::TrapezoidBounds>(
      minHalflengthX, maxHalflengthX, halflengthZ);

  // the Acts segmentation of the DigitizationModule
  size_t bins0 =
      (gridSizeX != 0)
          ? static_cast<size_t>((2 * maxHalflengthX) / (gridSizeX * scalor))
          : 0;
  size_t bins1 =
      (gridSizeZ != 0)
          ? static_cast<size_t>((2 * halflengthZ) / (gridSizeZ * scalor))
          : 0;

  std::shared_ptr<const Acts::CartesianSegmentation> actsSegmentation =
      std::make_shared<const Acts::CartesianSegmentation>(bounds, bins0, bins1);
  // finally create the digitization module
  // @todo set lorentz angle
  return (std::make_shared<const Acts::DigitizationModule>(actsSegmentation,
                                                           thickness, 1, 0));
}

dd4hep::xml::Component getNodeByStrAttr(const dd4hep::xml::Handle_t& mother,
                                        const std::string& nodeName,
                                        const std::string& attrName,
                                        const std::string& attrValue) {
  for (dd4hep::xml::Collection_t xCompColl(mother, nodeName.c_str());
       nullptr != xCompColl; ++xCompColl) {
    if (xCompColl.attr<std::string>(attrName.c_str()) == attrValue) {
      return static_cast<dd4hep::xml::Component>(xCompColl);
    }
  }
  // in case there was no xml daughter with matching name
  return dd4hep::xml::Component(nullptr);
}

double getAttrValueWithFallback(const dd4hep::xml::Component& node,
                                const std::string& attrName,
                                const double& defaultValue) {
  if (node.hasAttr(_Unicode(attrName.c_str()))) {
    return node.attr<double>(attrName.c_str());
  } else {
    return defaultValue;
  }
}
}  // namespace utils
}  // namespace det
