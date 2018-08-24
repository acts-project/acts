// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Units.hpp"
#include "DD4hep/CartesianGridXY.h"

Acts::DD4hepDetElement::DD4hepDetElement(
    const dd4hep::DetElement                     detElement,
    const std::string&                           axes,
    double                                       scalor,
    bool                                         isDisc,
    std::shared_ptr<const Acts::SurfaceMaterial> material,
    bool                                         buildDigitizationModules,
    std::shared_ptr<const DigitizationModule>    digiModule)
  : Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),
                              detElement.nominal().worldTransformation(),
                              detElement.placement().ptr(),
                              axes,
                              scalor,
                              isDisc,
                              material)
  , m_detElement(std::move(detElement))
  , m_digiModule(digiModule)

{
  // if wanted access the segmentation and create digitization module if not
  // handed over
  if (buildDigitizationModules && m_detElement.volume().isSensitive()
      && !m_digiModule) {
    dd4hep::SensitiveDetector sensDet(
        m_detElement.volume().sensitiveDetector());
    // assign the segmentation
    m_segmentation = sensDet.readout().segmentation();
    // @todo when sensitive detector is component

    // Now create the DigitizationModule
    // get the grid from DD4hep
    dd4hep::CartesianGridXY cartesianGrid = m_segmentation;
    if (cartesianGrid.isValid()) {
      // the segmentation of the DigitizationModule
      std::shared_ptr<const CartesianSegmentation> segmentation = nullptr;
      // find out the surface bounds
      // rectangular case
      std::shared_ptr<const RectangleBounds> bounds(
          dynamic_cast<const RectangleBounds*>(
              this->surface().bounds().clone()));
      if (bounds) {
        // the number of bins is the total length of the module divided by the
        // cell length
        size_t bins0 = (cartesianGrid.gridSizeX() != 0)
            ? (bounds->halflengthX() * 2.) / cartesianGrid.gridSizeX()
            : 0;
        size_t bins1 = (cartesianGrid.gridSizeY()) != 0.0
            ? (bounds->halflengthY() * 2.) / cartesianGrid.gridSizeY()
            : 0;
        // create the segmentation
        segmentation = std::make_shared<const CartesianSegmentation>(
            bounds, bins0, bins1);
      }
      // trapezoidal case
      std::shared_ptr<const TrapezoidBounds> tbounds(
          dynamic_cast<const TrapezoidBounds*>(
              this->surface().bounds().clone()));
      if (tbounds) {
        // the number of bins is the total length of the module divided by the
        // cell length
        size_t bins0 = (cartesianGrid.gridSizeX() != 0)
            ? (tbounds->maxHalflengthX() * 2.) / cartesianGrid.gridSizeX()
            : 0;
        size_t bins1 = (cartesianGrid.gridSizeY() != 0)
            ? (tbounds->halflengthY() * 2.) / cartesianGrid.gridSizeY()
            : 0;
        // create the segmentation
        segmentation = std::make_shared<const CartesianSegmentation>(
            tbounds, bins0, bins1);
      }
      // finally create the digitization module
      // @todo set lorentz angle
      m_digiModule = std::make_shared<const DigitizationModule>(
          segmentation, this->thickness(), 1, 0);
    }
  }
}

std::shared_ptr<const Acts::DigitizationModule>
Acts::DD4hepDetElement::digitizationModule() const
{
  return m_digiModule;
}
