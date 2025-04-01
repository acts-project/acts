// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoMaterialConverter.hpp"

#include "Acts/Material/Material.hpp"

#include "TGeoMaterial.h"

Acts::MaterialSlab Acts::TGeoMaterialConverter::materialSlab(
    const TGeoMaterial& tgMaterial, long double thicknessIn,
    long double thicknessOut, const Options& options) {
  // Scalable properties
  long double matX0 = tgMaterial.GetRadLen();
  long double matL0 = tgMaterial.GetIntLen();
  long double matRho = tgMaterial.GetDensity();

  // X0, L0, rho scale with the thickness
  long double cFactor = thicknessIn / thicknessOut;
  long double uScalor = options.unitLengthScalor;
  long double rScalar =
      options.unitMassScalor / pow(options.unitLengthScalor, 3);

  auto material = Material::fromMassDensity(
      matX0 * uScalor / cFactor, matL0 * uScalor / cFactor, tgMaterial.GetA(),
      tgMaterial.GetZ(), matRho * rScalar * cFactor);

  return MaterialSlab(material, thicknessOut);
}
