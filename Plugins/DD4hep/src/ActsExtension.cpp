// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include <boost/algorithm/string.hpp>
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Units.hpp"
#include "DD4hep/CartesianGridXY.h"

Acts::ActsExtension::ActsExtension(const Config& cfg)
  : Acts::IActsExtension(), m_material(nullptr)
{
  setConfiguration(cfg);
}

Acts::ActsExtension::ActsExtension(const ActsExtension& det,
                                   const dd4hep::DetElement& /*elem*/)
  : Acts::IActsExtension(), m_cfg(det.m_cfg), m_material(det.m_material)
{
}

Acts::ActsExtension::ActsExtension(
    std::shared_ptr<const DigitizationModule> digiModule)
  : Acts::IActsExtension()
  , m_material(nullptr)
  , m_digitizationModule(std::move(digiModule))
{
}

Acts::ActsExtension::ActsExtension(
    const std::vector<std::pair<dd4hep::Material, double>>& materials,
    std::shared_ptr<const DigitizationModule> digiModule)
  : Acts::IActsExtension()
  , m_material(nullptr)
  , m_digitizationModule(std::move(digiModule))
{
  std::vector<Acts::MaterialProperties> partialMaterial;
  partialMaterial.reserve(materials.size());
  for (auto& mat : materials) {
    Acts::Material pm{float(mat.first.radLength() * units::_cm),
                      float(mat.first.intLength() * units::_cm),
                      float(mat.first.A()),
                      float(mat.first.Z()),
                      float(mat.first.density() / pow(Acts::units::_cm, 3))};
    partialMaterial.push_back(
        Acts::MaterialProperties(pm, mat.second * units::_mm));
  }
  //  Create homogenous surface material with averaged material properties
  Acts::MaterialProperties matprop(partialMaterial);
  m_material = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matprop);
}

void
Acts::ActsExtension::setConfiguration(const Acts::ActsExtension::Config& config)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = config;
}
