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

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_CASE(AccumulatedMaterialMaterial_test)
  {
    AccumulatedVolumeMaterial avm;

    // Try averaging without content
    Material result = avm.average();
    BOOST_CHECK_EQUAL(result, Material());

    // Test that the mean of a single material is the material itself
    Material mat1(1., 2., 3., 4., 5.);
    avm.accumulate(mat1);
    result = avm.average();
    BOOST_CHECK_EQUAL(result.X0(), mat1.X0());
    BOOST_CHECK_EQUAL(result.L0(), mat1.L0());
    BOOST_CHECK_EQUAL(result.A(), mat1.A());
    BOOST_CHECK_EQUAL(result.Z(), mat1.Z());
    BOOST_CHECK_EQUAL(result.rho(), mat1.rho());

    // Test that the mean is actually calculated
    Material mat2(6., 7., 8., 9., 10.);
    avm.accumulate(mat2);
    result = avm.average();
    BOOST_CHECK_EQUAL(result.X0(), 0.5 * (mat1.X0() + mat2.X0()));
    BOOST_CHECK_EQUAL(result.L0(), 0.5 * (mat1.L0() + mat2.L0()));
    BOOST_CHECK_EQUAL(result.A(), 0.5 * (mat1.A() + mat2.A()));
    BOOST_CHECK_EQUAL(result.Z(), 0.5 * (mat1.Z() + mat2.Z()));
    BOOST_CHECK_EQUAL(result.rho(), 0.5 * (mat1.rho() + mat2.rho()));

    // Test that the mean of vacuums is a vacuum
    AccumulatedVolumeMaterial avm2;
    avm2.accumulate(Material());
    avm2.accumulate(Material());
    result = avm2.average();
    BOOST_CHECK_EQUAL(result, Material());

    // Add vacuum to the material and test it
    avm.accumulate(Material());
    result = avm.average();
    BOOST_CHECK_EQUAL(result.X0(), 0.25 * (mat1.X0() + mat2.X0()) * 3);
    BOOST_CHECK_EQUAL(result.L0(), 0.25 * (mat1.L0() + mat2.L0()) * 3);
    BOOST_CHECK_EQUAL(result.A(), (mat1.A() + mat2.A()) / 3);
    BOOST_CHECK_EQUAL(result.Z(), (mat1.Z() + mat2.Z()) / 3);
    BOOST_CHECK_EQUAL(result.rho(), (mat1.rho() + mat2.rho()) / 3);

    // Add material to the vacuum and test it
    avm2.accumulate(mat1);
    result = avm2.average();
    BOOST_CHECK_EQUAL(result.X0(), mat1.X0() * 3);
    BOOST_CHECK_EQUAL(result.L0(), mat1.L0() * 3);
    BOOST_CHECK_EQUAL(result.A(), mat1.A() / 3);
    BOOST_CHECK_EQUAL(result.Z(), mat1.Z() / 3);
    BOOST_CHECK_EQUAL(result.rho(), mat1.rho() / 3);
  }
}  // namespace Test
}  // namespace Acts