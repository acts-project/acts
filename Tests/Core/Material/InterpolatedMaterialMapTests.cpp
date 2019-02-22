// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedBFieldMap_tests.cpp

// clang-format off
#define BOOST_TEST_MODULE Mapped material tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include <array>

namespace Acts {

namespace Test {

  constexpr unsigned int dim = 2;

  ActsVectorD<dim>
  trafoGlobalToLocal(const Vector3D& global)
  {
    return {global.x(), global.y()};
  }

  BOOST_AUTO_TEST_CASE(InterpolatedMaterialMap_test)
  {
    std::array<double, dim> lowerLeft{{0., 0.}};
    std::array<double, dim> upperRight{{1., 1.}};
    Material mat(1., 2., 3., 4., 5.);

    InterpolatedMaterialMap::MaterialCell<dim> materialCell(
        trafoGlobalToLocal, lowerLeft, upperRight, mat);
  }
}  // namespace Test

}  // namespace Acts
