// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoMaterialConverter.hpp"

#include "Acts/Material/Material.hpp"

#include "TGeoMaterial.h"

Acts::MaterialSlab Acts::TGeoMaterialConverter::materialSlab(
    const TGeoMaterial& tgMaterial, ActsScalar thicknessIn,
    ActsScalar thicknessOut, const Options& options) {
  ActsScalar matA = tgMaterial.GetA();
  ActsScalar matZ = tgMaterial.GetZ();
  ActsScalar matX0 = tgMaterial.GetRadLen();
  ActsScalar matL0 = tgMaterial.GetIntLen();
  // @todo TGeo to Acts density conversion missing
  ActsScalar matRho = tgMaterial.GetDensity();

  // X0, L0, rho scale with the thicknes
  ActsScalar cFactor = thicknessIn / thicknessOut;
  ActsScalar uScalor = options.unitLengthScalor;
  ActsScalar rScalar =
      options.unitMassScalor / pow(options.unitLengthScalor, 3);
  auto material = Material::fromMassDensity(matX0 * uScalor / cFactor,
                                            matL0 * uScalor / cFactor, matA,
                                            matZ, matRho * rScalar * cFactor);

  return MaterialSlab(std::move(material), thicknessOut);
}
