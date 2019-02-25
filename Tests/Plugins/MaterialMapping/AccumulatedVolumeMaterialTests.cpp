// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AccumulatedVolumeMaterial Tests
#include <boost/test/included/unit_test.hpp>
#include <vector>
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedVolumeMaterial.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_CASE(AccumulatedMaterialMaterial_test)
  {
    AccumulatedVolumeMaterial avm;

    // Try averaging without content
    Material result = avm.average();
    CHECK_CLOSE_REL(result, Material(), 1e-4);

    // Test that the mean of a single material is the material itself
    Material mat1(1., 2., 3., 4., 5.);
    avm.accumulate(mat1);
    result = avm.average();
    CHECK_CLOSE_REL(result.X0(), mat1.X0(), 1e-4);
    CHECK_CLOSE_REL(result.L0(), mat1.L0(), 1e-4);
    CHECK_CLOSE_REL(result.A(), mat1.A(), 1e-4);
    CHECK_CLOSE_REL(result.Z(), mat1.Z(), 1e-4);
    CHECK_CLOSE_REL(result.rho(), mat1.rho(), 1e-4);

    // Test that the mean is actually calculated
    Material mat2(6., 7., 8., 9., 10.);
    avm.accumulate(mat2);
    result = avm.average();
    CHECK_CLOSE_REL(result.X0(), 0.5 * (mat1.X0() + mat2.X0()), 1e-4);
    CHECK_CLOSE_REL(result.L0(), 0.5 * (mat1.L0() + mat2.L0()), 1e-4);
    CHECK_CLOSE_REL(result.A(), 0.5 * (mat1.A() + mat2.A()), 1e-4);
    CHECK_CLOSE_REL(result.Z(), 0.5 * (mat1.Z() + mat2.Z()), 1e-4);
    CHECK_CLOSE_REL(result.rho(), 0.5 * (mat1.rho() + mat2.rho()), 1e-4);

    // Test that the mean of vacuums is a vacuum
    AccumulatedVolumeMaterial avm2;
    avm2.accumulate(Material());
    avm2.accumulate(Material());
    result = avm2.average();
    CHECK_CLOSE_REL(result, Material(), 1e-4);

    // Add vacuum to the material and test it
    avm.accumulate(Material());
    result = avm.average();
    CHECK_CLOSE_REL(result.X0(), 0.25 * (mat1.X0() + mat2.X0()) * 3, 1e-4);
    CHECK_CLOSE_REL(result.L0(), 0.25 * (mat1.L0() + mat2.L0()) * 3, 1e-4);
    CHECK_CLOSE_REL(result.A(), (mat1.A() + mat2.A()) / 3, 1e-4);
    CHECK_CLOSE_REL(result.Z(), (mat1.Z() + mat2.Z()) / 3, 1e-4);
    CHECK_CLOSE_REL(result.rho(), (mat1.rho() + mat2.rho()) / 3, 1e-4);

    // Add material to the vacuum and test it
    avm2.accumulate(mat1);
    result = avm2.average();
    CHECK_CLOSE_REL(result.X0(), mat1.X0() * 3, 1e-4);
    CHECK_CLOSE_REL(result.L0(), mat1.L0() * 3, 1e-4);
    CHECK_CLOSE_REL(result.A(), mat1.A() / 3, 1e-4);
    CHECK_CLOSE_REL(result.Z(), mat1.Z() / 3, 1e-4);
    CHECK_CLOSE_REL(result.rho(), mat1.rho() / 3, 1e-4);
  }
}  // namespace Test
}  // namespace Acts