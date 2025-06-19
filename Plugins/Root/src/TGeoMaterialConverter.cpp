// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/TGeoMaterialConverter.hpp"

#include "Acts/Material/Material.hpp"

#include "TGeoMaterial.h"

Acts::MaterialSlab Acts::TGeoMaterialConverter::materialSlab(
    const TGeoMaterial& tgMaterial, double thicknessIn, double thicknessOut,
    const Options& options) {
  // Scalable properties
  double matX0 = tgMaterial.GetRadLen();
  double matL0 = tgMaterial.GetIntLen();
  double matRho = tgMaterial.GetDensity();

  // X0, L0, rho scale with the thickness
  double cFactor = thicknessIn / thicknessOut;
  double uScalor = options.unitLengthScalor;
  double rScalar = options.unitMassScalor / pow(options.unitLengthScalor, 3);

  auto material =
      Material::fromMassDensity(static_cast<float>(matX0 * uScalor / cFactor),
                                static_cast<float>(matL0 * uScalor / cFactor),
                                static_cast<float>(tgMaterial.GetA()),
                                static_cast<float>(tgMaterial.GetZ()),
                                static_cast<float>(matRho * rScalar * cFactor));

  return MaterialSlab(material, static_cast<float>(thicknessOut));
}
