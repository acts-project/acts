// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hepPlugins/ActsExtension.hpp"
#include <boost/algorithm/string.hpp>
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Units.hpp"
#include "DD4hep/CartesianGridXY.h"

std::shared_ptr<const Acts::DigitizationModule>
Acts::rectangleDigiModule(double                      halflengthX,
                          double                      halflengthY,
                          double                      thickness,
                          const dd4hep::Segmentation& segmentation)
{
  // convert to Acts units
  double scalor = units::_cm;
  halflengthX *= scalor;
  halflengthY *= scalor;
  thickness *= scalor;

  auto bounds
      = std::make_shared<const RectangleBounds>(halflengthX, halflengthY);
  dd4hep::CartesianGridXY cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
        ? halflengthX / cartesianGrid.gridSizeX()
        : 0;
    size_t bins1 = (cartesianGrid.gridSizeY() != 0)
        ? halflengthY / cartesianGrid.gridSizeY()
        : 0;

    std::shared_ptr<const CartesianSegmentation> actsSegmentation
        = std::make_shared<const CartesianSegmentation>(bounds, bins0, bins1);
    // finally create the digitization module
    // @todo set lorentz angle
    return (std::make_shared<const DigitizationModule>(
        actsSegmentation, thickness, 1, 0));
  }
  return nullptr;
}

std::shared_ptr<const Acts::DigitizationModule>
Acts::trapezoidalDigiModule(double                      minHalflengthX,
                            double                      maxHalflengthX,
                            double                      halflengthY,
                            double                      thickness,
                            const dd4hep::Segmentation& segmentation)
{
  // convert to Acts units
  double scalor = units::_cm;
  minHalflengthX *= scalor;
  maxHalflengthX *= scalor;
  halflengthY *= scalor;
  thickness *= scalor;

  auto bounds = std::make_shared<const TrapezoidBounds>(
      minHalflengthX, maxHalflengthX, halflengthY);
  ;
  dd4hep::CartesianGridXY cartesianGrid = segmentation;
  if (cartesianGrid.isValid()) {
    // the Acts segmentation of the DigitizationModule
    size_t bins0 = (cartesianGrid.gridSizeX() != 0)
        ? maxHalflengthX / cartesianGrid.gridSizeX()
        : 0;
    size_t bins1 = (cartesianGrid.gridSizeY() != 0)
        ? halflengthY / cartesianGrid.gridSizeY()
        : 0;

    std::shared_ptr<const CartesianSegmentation> actsSegmentation
        = std::make_shared<const CartesianSegmentation>(bounds, bins0, bins1);
    // finally create the digitization module
    // @todo set lorentz angle
    return (std::make_shared<const DigitizationModule>(
        actsSegmentation, thickness, 1, 0));
  }
  return nullptr;
}

Acts::ActsExtension::ActsExtension(const Config& cfg)
  : Acts::IActsExtension(), m_material(nullptr), m_digiModule(nullptr)
{
  setConfiguration(cfg);
}

Acts::ActsExtension::ActsExtension(
    const std::vector<std::pair<dd4hep::Material, double>>& materials,
    std::shared_ptr<const DigitizationModule> digiModule)
  : Acts::IActsExtension(), m_material(nullptr), m_digiModule(digiModule)
{
  Acts::MaterialProperties matprop;
  for (auto& mat : materials) {
    matprop.add(
        MaterialProperties(mat.first.radLength() * units::_cm,
                           mat.first.intLength() * units::_cm,
                           mat.first.A(),
                           mat.first.Z(),
                           mat.first.density() / pow(Acts::units::_cm, 3),
                           mat.second * units::_mm));
  }

  //  Create homogenous surface material with averaged material properties
  m_material = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matprop);
}

Acts::ActsExtension::ActsExtension(
    std::shared_ptr<const DigitizationModule> digiModule)
  : Acts::IActsExtension(), m_material(nullptr), m_digiModule(digiModule)
{
}

Acts::ActsExtension::ActsExtension(const ActsExtension& det,
                                   const dd4hep::DetElement&)
  : Acts::IActsExtension()
  , m_cfg(det.m_cfg)
  , m_material(det.m_material)
  , m_digiModule(det.m_digiModule)
{
}

void
Acts::ActsExtension::setConfiguration(const Acts::ActsExtension::Config& config)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = config;
}
