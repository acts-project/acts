// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE ProtoLayerTests
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>

#include <cmath>
#include "Acts/Layers/ProtoLayer.hpp"

namespace Acts {

namespace Test {
  namespace Layers {
    BOOST_AUTO_TEST_SUITE(Layers)

    BOOST_AUTO_TEST_CASE(ProtoLayer_radialDistance)
    {
      ProtoLayer pl;

      Vector3D p1(10, 0, 0);
      Vector3D p2(0, 10, 0);
      BOOST_CHECK_CLOSE_FRACTION(pl.radialDistance(p1, p2), sqrt(50), 1e-6);

      Vector3D p3(-5, 5, 0);
      Vector3D p4(5, 5, 0);
      BOOST_CHECK_CLOSE_FRACTION(pl.radialDistance(p3, p4), 5, 1e-6);

      Vector3D p5(6, 6, 0);
      Vector3D p6(8, 9, 0);
      BOOST_CHECK_CLOSE_FRACTION(pl.radialDistance(p5, p6),
                                 sqrt(6 * 6 + 6 * 6),
                                 1e-6);

      Vector3D p7(0, 10, 0);
      Vector3D p8(5, 5, 0);
      BOOST_CHECK_CLOSE_FRACTION(pl.radialDistance(p7, p8), sqrt(50), 1e-6);

      Vector3D p9(13, 2, 0);
      Vector3D p10(13, -2, 0);
      BOOST_CHECK_CLOSE_FRACTION(pl.radialDistance(p9, p10), 13, 1e-6);
    }

    BOOST_AUTO_TEST_SUITE_END()
  }  // namespace Layers
}  // namespace Test

}  // namespace Acts
