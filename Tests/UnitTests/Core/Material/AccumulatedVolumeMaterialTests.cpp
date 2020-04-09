// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/AccumulatedVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(accumulated_material)

BOOST_AUTO_TEST_CASE(vacuum) {
  AccumulatedVolumeMaterial avm;

  // averaging over nothing is vacuum
  BOOST_TEST(avm.average() == Material());

  // averaging over vacuum is still vacuum
  avm.accumulate(MaterialProperties(1));
  BOOST_TEST(avm.average() == Material());
}

BOOST_AUTO_TEST_CASE(single_material) {
  Material mat(1., 2., 3., 4., 5.);
  MaterialProperties matprop(mat, 1);
  AccumulatedVolumeMaterial avm;
  // mean of a single material should be the same material again for a thickness
  // of 1
  avm.accumulate(matprop);
  {
    auto result = avm.average();
    CHECK_CLOSE_REL(result.X0(), mat.X0(), 1e-4);
    CHECK_CLOSE_REL(result.L0(), mat.L0(), 1e-4);
    CHECK_CLOSE_REL(result.Ar(), mat.Ar(), 1e-4);
    CHECK_CLOSE_REL(result.Z(), mat.Z(), 1e-4);
    CHECK_CLOSE_REL(result.massDensity(), mat.massDensity(), 1e-4);
  }
  // adding a vacuum step changes the average
  avm.accumulate(MaterialProperties(1));
  {
    auto result = avm.average();
    // less scattering in vacuum, larger radiation length
    CHECK_CLOSE_REL(result.X0(), 2 * mat.X0(), 1e-4);
    CHECK_CLOSE_REL(result.L0(), 2 * mat.L0(), 1e-4);
    // less material, lower density
    CHECK_CLOSE_REL(result.Ar(), 0.5 * mat.Ar(), 1e-4);
    CHECK_CLOSE_REL(result.Z(), 0.5 * mat.Z(), 1e-4);
    CHECK_CLOSE_REL(result.massDensity(), 0.5 * mat.massDensity(), 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(two_materials) {
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);

  MaterialProperties matprop1(mat1, 1);
  MaterialProperties matprop2(mat2, 1);

  AccumulatedVolumeMaterial avm;
  avm.accumulate(matprop1);
  avm.accumulate(matprop2);
  auto result = avm.average();
  CHECK_CLOSE_REL(result.X0(), 2. / (1. / 1. + 1. / 6.), 1e-4);
  CHECK_CLOSE_REL(result.L0(), 2. / (1. / 2. + 1. / 7.), 1e-4);
  CHECK_CLOSE_REL(result.Ar(), 0.5 * (3. + 8.), 1e-4);
  CHECK_CLOSE_REL(result.Z(), 0.5 * (4. + 9.), 1e-4);
  CHECK_CLOSE_REL(result.massDensity(), 0.5 * (5. + 10.), 1e-4);
}

BOOST_AUTO_TEST_CASE(two_materials_different_lengh) {
  Material mat1(1., 2., 3., 4., 5.);
  Material mat2(6., 7., 8., 9., 10.);

  MaterialProperties matprop1(mat1, 0.5);
  MaterialProperties matprop2(mat2, 2);

  AccumulatedVolumeMaterial avm;
  avm.accumulate(matprop1);
  avm.accumulate(matprop2);
  auto result = avm.average();
  CHECK_CLOSE_REL(result.X0(), 2.5 / (0.5 / 1. + 2. / 6.), 1e-4);
  CHECK_CLOSE_REL(result.L0(), 2.5 / (0.5 / 2. + 2. / 7.), 1e-4);
  CHECK_CLOSE_REL(result.Ar(), 0.5 * (3. + 8.), 1e-4);
  CHECK_CLOSE_REL(result.Z(), 0.5 * (4. + 9.), 1e-4);
  CHECK_CLOSE_REL(result.massDensity(), 0.5 * (5. + 10.), 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
