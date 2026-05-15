// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"

namespace ActsTests {

inline Acts::Material makeBeryllium() {
  using namespace Acts::UnitLiterals;
  return Acts::Material::fromMolarDensity(35.28_cm, 42.10_cm, 9.012, 4,
                                          (1.848 / 9.012) * 1_mol / 1_cm3);
}

inline Acts::Material makeSilicon() {
  using namespace Acts::UnitLiterals;
  return Acts::Material::fromMolarDensity(9.370_cm, 46.52_cm, 28.0855, 14,
                                          (2.329 / 28.0855) * 1_mol / 1_cm3);
}

inline Acts::Material makeLiquidArgon() {
  using namespace Acts::UnitLiterals;
  return Acts::Material::fromMassDensity(140.0343868497_mm, 857.0639538668_mm,
                                         39.9476933511, 18,
                                         1.396 * 1_g / 1_cm3);
}

inline Acts::Material makeIron() {
  using namespace Acts::UnitLiterals;
  return Acts::Material::fromMassDensity(17.57493465097_mm, 169.9030027586_mm,
                                         55.845110798, 26, 7.874 * 1_g / 1_cm3);
}

/// Build material slab corresponding to 1 radiation and interaction length.
inline Acts::MaterialSlab makeUnitSlab() {
  using namespace Acts::UnitLiterals;
  // silicon-like material with higher X0 and lower L0
  return {Acts::Material::fromMolarDensity(20_cm, 20_cm, 28.0855, 14,
                                           (2.329 / 28.0855) * 1_mol / 1_cm3),
          20_cm};
}

/// Build material slab corresponding to 1% of radiation and interaction length.
inline Acts::MaterialSlab makePercentSlab() {
  auto slab = makeUnitSlab();
  slab.scaleThickness(0.01);
  return slab;
}

}  // namespace ActsTests
