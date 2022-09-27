// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

namespace Acts {
namespace Experimental {
class DetectorVolume {};
}  // namespace Experimental
}  // namespace Acts

using namespace Acts::Experimental;

// A test context
Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

DetectorVolume dVolume;
auto dTransform = Acts::Transform3::Identity();
auto pGenerator = detail::defaultPortalGenerator();

BOOST_AUTO_TEST_CASE(CylindricalPortalGenerator) {
  // Access Vectors, they should yield the detector volume
  // as a return on nextVolume(...) call
  Acts::Vector3 negPos(-200., 0., 0.);
  Acts::Vector3 negDir(0., 0., 1.);

  Acts::Vector3 posPos(200., 0., 0.);
  Acts::Vector3 posDir(0., 0., -1.);

  Acts::Vector3 outerPos(100., 0., 0.);
  Acts::Vector3 outerDir(-1., 0., 0.);

  Acts::Vector3 innerPos(10., 0., 0.);
  Acts::Vector3 innerDir(1., 0., 0.);

  // Filled Cylinder
  Acts::CylinderVolumeBounds cBar(0., 100, 200.);

  auto cBarPortals = detail::portals(dTransform, cBar, dVolume);
  
  BOOST_CHECK(cBarPortals.size() == 3u);
  // Check they are not nullptrs
  for (const auto& p : cBarPortals) {
    BOOST_CHECK(p != nullptr);
  }

  // Pointing inside the volume
  BOOST_CHECK(cBarPortals[0]->nextVolume(tContext, negPos, negDir) == &dVolume);
  BOOST_CHECK(cBarPortals[1]->nextVolume(tContext, posPos, posDir) == &dVolume);
  BOOST_CHECK(cBarPortals[2]->nextVolume(tContext, outerPos, outerDir) ==
              &dVolume);
  // Pointing to nowhere land
  BOOST_CHECK(cBarPortals[0]->nextVolume(tContext, negPos, -negDir) == nullptr);
  BOOST_CHECK(cBarPortals[1]->nextVolume(tContext, posPos, -posDir) == nullptr);
  BOOST_CHECK(cBarPortals[2]->nextVolume(tContext, outerPos, -outerDir) ==
              nullptr);

  // Tube Cylinder
  Acts::CylinderVolumeBounds cTube(10., 100, 200.);
  auto cTubePortals = detail::portals(dTransform, cTube, dVolume);
  BOOST_CHECK(cTubePortals.size() == 4u);
  // Check they are not nullptrs
  for (const auto& p : cTubePortals) {
    BOOST_CHECK(p != nullptr);
  }

  // Pointing inside the volume
  BOOST_CHECK(cTubePortals[0]->nextVolume(tContext, negPos, negDir) ==
              &dVolume);
  BOOST_CHECK(cTubePortals[1]->nextVolume(tContext, posPos, posDir) ==
              &dVolume);
  BOOST_CHECK(cTubePortals[2]->nextVolume(tContext, outerPos, outerDir) ==
              &dVolume);
  BOOST_CHECK(cTubePortals[3]->nextVolume(tContext, innerPos, innerDir) ==
              &dVolume);

  // Pointing to nowhere land
  BOOST_CHECK(cTubePortals[0]->nextVolume(tContext, negPos, -negDir) ==
              nullptr);
  BOOST_CHECK(cTubePortals[1]->nextVolume(tContext, posPos, -posDir) ==
              nullptr);
  BOOST_CHECK(cTubePortals[2]->nextVolume(tContext, outerPos, -outerDir) ==
              nullptr);
  BOOST_CHECK(cTubePortals[3]->nextVolume(tContext, innerPos, -innerDir) ==
              nullptr);

  // Sectoral tube cylinder
  Acts::ActsScalar alpha = 0.25 * M_PI;
  Acts::ActsScalar r = 50;

  Acts::Vector3 negPhiSecPos(r * std::cos(-alpha), r * std::sin(-alpha), 0.);
  Acts::Vector3 negPhiSecDir(-r * std::cos(-alpha), r * std::sin(-alpha), 0.);
  Acts::Vector3 posPhiSecPos(r * std::cos(alpha), r * std::sin(alpha), 0.);
  Acts::Vector3 posPhiSecDir(r * std::cos(alpha), -r * std::sin(alpha), 0.);

  Acts::CylinderVolumeBounds cTubeSector(10., 100., 200., alpha, 0.);
  auto cTubeSectorPortals = detail::portals(dTransform, cTubeSector, dVolume);
  BOOST_CHECK(cTubeSectorPortals.size() == 6u);
  // Check they are not nullptrs
  for (const auto& p : cTubeSectorPortals) {
    BOOST_CHECK(p != nullptr);
  }
  // Pointing inside the volume
  BOOST_CHECK(cTubeSectorPortals[0]->nextVolume(tContext, negPos, negDir) ==
              &dVolume);
  BOOST_CHECK(cTubeSectorPortals[1]->nextVolume(tContext, posPos, posDir) ==
              &dVolume);
  BOOST_CHECK(cTubeSectorPortals[2]->nextVolume(tContext, outerPos, outerDir) ==
              &dVolume);
  BOOST_CHECK(cTubeSectorPortals[3]->nextVolume(tContext, innerPos, innerDir) ==
              &dVolume);
  BOOST_CHECK(cTubeSectorPortals[4]->nextVolume(tContext, negPhiSecPos,
                                                negPhiSecDir) == &dVolume);
  BOOST_CHECK(cTubeSectorPortals[5]->nextVolume(tContext, posPhiSecPos,
                                                posPhiSecDir) == &dVolume);

  // Pointing to nowhere land
  BOOST_CHECK(cTubeSectorPortals[0]->nextVolume(tContext, negPos, -negDir) ==
              nullptr);
  BOOST_CHECK(cTubeSectorPortals[1]->nextVolume(tContext, posPos, -posDir) ==
              nullptr);
  BOOST_CHECK(cTubeSectorPortals[2]->nextVolume(tContext, outerPos,
                                                -outerDir) == nullptr);
  BOOST_CHECK(cTubeSectorPortals[3]->nextVolume(tContext, innerPos,
                                                -innerDir) == nullptr);
  BOOST_CHECK(cTubeSectorPortals[4]->nextVolume(tContext, negPhiSecPos,
                                                -negPhiSecDir) == nullptr);
  BOOST_CHECK(cTubeSectorPortals[5]->nextVolume(tContext, posPhiSecPos,
                                                -posPhiSecDir) == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
