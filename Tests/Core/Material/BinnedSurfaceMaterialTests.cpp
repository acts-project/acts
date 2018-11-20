// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE SurfaceMaterial Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

namespace Test {

  /// Test the constructors
  BOOST_AUTO_TEST_CASE(BinnedSurfaceMaterial_construction_test)
  {

    BinUtility xyBinning(2, -1., 1., open, binX);
    xyBinning += BinUtility(3, -3., 3., open, binY);

    // Constructor a few material properties
    MaterialProperties a00(1., 2., 3., 4., 5., 6.);
    MaterialProperties a01(2., 3., 4., 5., 6., 7.);
    MaterialProperties a02(3., 4., 5., 6., 7., 8.);
    MaterialProperties a10(4., 5., 6., 7., 8., 9.);
    MaterialProperties a11(5., 6., 7., 8., 9., 10.);
    MaterialProperties a12(6., 7., 8., 9., 10., 11.);

    // Prepare the matrix
    std::vector<MaterialProperties> l0 = {std::move(a00), std::move(a10)};
    std::vector<MaterialProperties> l1 = {std::move(a01), std::move(a11)};
    std::vector<MaterialProperties> l2 = {std::move(a02), std::move(a12)};

    // Build the matrix
    std::vector<std::vector<MaterialProperties>> m
        = {std::move(l0), std::move(l1), std::move(l2)};

    // Create the material
    BinnedSurfaceMaterial bsm(xyBinning, std::move(m));

    // Copy the material
    BinnedSurfaceMaterial bsmCopy(bsm);

    // Assignment operation
    BinnedSurfaceMaterial bsmAssigned = bsm;

    // Move constructor
    BinnedSurfaceMaterial bsmMoved(std::move(bsmCopy));

    // Move assigned
    BinnedSurfaceMaterial bsmMoveAssigned(std::move(bsmAssigned));
  }

}  // namespace Test
}  // namespace Acts
