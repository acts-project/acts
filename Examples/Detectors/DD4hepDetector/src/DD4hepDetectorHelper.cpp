// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetectorHelper.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cstddef>

#include "DD4hep/CartesianGridXY.h"
#include "DD4hep/Segmentations.h"

using namespace ActsExamples::DD4hep;

std::shared_ptr<const Acts::DigitizationModule>
DD4hepDetectorHelper::rectangleDigiModule(
    double halflengthX, double halflengthY, double thickness,
    const dd4hep::Segmentation& segmentation) {
  // convert to Acts units
  double scalor = Acts::UnitConstants::cm;
  halflengthX *= scalor;
  halflengthY *= scalor;
  thickness *= scalor;

  auto bounds =
      std::make_shared<const Acts::RectangleBounds>(halflengthX, halflengthY);
  dd4hep::CartesianGridXY cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    double gridSizeX = cartesianGrid.gridSizeX() * scalor;
    double gridSizeY = cartesianGrid.gridSizeY() * scalor;
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
                       ? static_cast<size_t>((2 * halflengthX) / gridSizeX)
                       : 0;
    size_t bins1 = (cartesianGrid.gridSizeY() != 0)
                       ? static_cast<size_t>((2 * halflengthY) / gridSizeY)
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

std::shared_ptr<const Acts::DigitizationModule>
DD4hepDetectorHelper::trapezoidalDigiModule(
    double minHalflengthX, double maxHalflengthX, double halflengthY,
    double thickness, const dd4hep::Segmentation& segmentation) {
  // convert to Acts units
  double scalor = Acts::UnitConstants::cm;
  minHalflengthX *= scalor;
  maxHalflengthX *= scalor;
  halflengthY *= scalor;
  thickness *= scalor;

  auto bounds = std::make_shared<const Acts::TrapezoidBounds>(
      minHalflengthX, maxHalflengthX, halflengthY);
  ;
  dd4hep::CartesianGridXY cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    double gridSizeX = cartesianGrid.gridSizeX() * scalor;
    double gridSizeY = cartesianGrid.gridSizeY() * scalor;
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
                       ? static_cast<size_t>((2 * maxHalflengthX) / gridSizeX)
                       : 0;
    size_t bins1 = (cartesianGrid.gridSizeY() != 0)
                       ? static_cast<size_t>((2 * halflengthY) / gridSizeY)
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
