// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AccumulatedVolumeMaterial Tests
#define BOOST_TEST_DYN_LINK
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
  avm.accumulate(Material());
  avm.accumulate(Material());
  BOOST_TEST(avm.average() == Material());
}

BOOST_AUTO_TEST_CASE(single_material) {
  Material mat(1., 2., 3., 4., 5.);

  AccumulatedVolumeMaterial avm;
  // mean of a single material should be the same material again
  avm.accumulate(mat);
  {
    auto result = avm.average();
    CHECK_CLOSE_REL(result.X0(), mat.X0(), 1e-4);
    CHECK_CLOSE_REL(result.L0(), mat.L0(), 1e-4);
    CHECK_CLOSE_REL(result.Ar(), mat.Ar(), 1e-4);
    CHECK_CLOSE_REL(result.Z(), mat.Z(), 1e-4);
    CHECK_CLOSE_REL(result.massDensity(), mat.massDensity(), 1e-4);
  }
  // adding a vacuum step changes the average
  avm.accumulate(Material());
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

  AccumulatedVolumeMaterial avm;
  avm.accumulate(mat1);
  avm.accumulate(mat2);
  auto result = avm.average();
  CHECK_CLOSE_REL(result.X0(), 0.5 * (1. + 6.), 1e-4);
  CHECK_CLOSE_REL(result.L0(), 0.5 * (2. + 7.), 1e-4);
  CHECK_CLOSE_REL(result.Ar(), 0.5 * (3. + 8.), 1e-4);
  CHECK_CLOSE_REL(result.Z(), 0.5 * (4. + 9.), 1e-4);
  CHECK_CLOSE_REL(result.massDensity(), 0.5 * (5. + 10.), 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
