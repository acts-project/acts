// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/context.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

using namespace Acts::UnitLiterals;

namespace Acts::Test {
GeometryContext gctx;

BOOST_AUTO_TEST_SUITE(PortalShellTests)

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  // - Cylinder
  // |           | no phi | phi |
  // | --------- | ------ | --- |
  // | rMin > 0  | 1      | 3   |
  // | rMin == 0 | 2      | 4   |

  std::size_t i = 1;
  auto makeVolume = [&](auto&&... pars) {
    TrackingVolume vol(Transform3::Identity(),
                       std::make_shared<CylinderVolumeBounds>(
                           std::forward<decltype(pars)>(pars)...));
    vol.setVolumeName("cyl" + std::to_string(i++));
    return std::move(vol);
  };

  auto cyl1 = makeVolume(30_mm, 40_mm, 100_mm);
  auto cyl2 = makeVolume(0_mm, 40_mm, 100_mm);
  auto cyl3 = makeVolume(30_mm, 40_mm, 100_mm, 45_degree);
  auto cyl4 = makeVolume(0_mm, 40_mm, 100_mm, 45_degree);

  SingleCylinderPortalShell shell1{cyl1};
  BOOST_CHECK_EQUAL(shell1.size(), 4);

  using enum CylinderPortalShell::Face;

  const auto* pDisc = shell1.portal(PositiveDisc);
  BOOST_REQUIRE_NE(pDisc, nullptr);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         -Vector3::UnitZ()),
                    &cyl1);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         Vector3::UnitZ()),
                    nullptr);

  const auto* nDisc = shell1.portal(NegativeDisc);
  BOOST_REQUIRE_NE(nDisc, nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         Vector3::UnitZ()),
                    &cyl1);

  SingleCylinderPortalShell shell2{cyl2};
  BOOST_CHECK_EQUAL(shell2.size(), 3);
  SingleCylinderPortalShell shell3{cyl3};
  BOOST_CHECK_EQUAL(shell3.size(), 6);
  SingleCylinderPortalShell shell4{cyl4};
  BOOST_CHECK_EQUAL(shell4.size(), 5);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
