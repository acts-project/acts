// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <exception>
#include <memory>

using namespace Acts::Experimental;

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(CylindricalDetectorVolume) {
  Acts::ActsScalar rInner = 10.;
  Acts::ActsScalar rOuter = 100.;
  Acts::ActsScalar zHalfL = 200.;

  Acts::Transform3 nominal = Acts::Transform3::Identity();

  auto fullCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(0., rOuter, zHalfL);

  auto portalGenerator = detail::defaultPortalGenerator();

  // Misconfigured - null pointer
  BOOST_CHECK_THROW(
      DetectorVolume::makeShared(tContext, nominal, nullptr, portalGenerator,
                                 "MisconfiguredFullCylinderVolume"),
      std::invalid_argument);

  // A full cylinder
  auto fullCylinderVolume = DetectorVolume::makeShared(
      tContext, nominal, std::move(fullCylinderBounds), portalGenerator,
      "FullCylinderVolume");

  BOOST_CHECK(fullCylinderVolume->surfaces().size() == 0u);
  BOOST_CHECK(fullCylinderVolume->volumes().size() == 0u);
  BOOST_CHECK(fullCylinderVolume->portals().size() == 3u);

  // A tube cylinder
  auto tubeCylinderBounds =
      std::make_unique<Acts::CylinderVolumeBounds>(rInner, rOuter, zHalfL);

  auto tubeCylinderVolume = DetectorVolume::makeShared(
      tContext, nominal, std::move(tubeCylinderBounds), portalGenerator,
      "TubeCylinderVolume");

  BOOST_CHECK(tubeCylinderVolume->surfaces().size() == 0u);
  BOOST_CHECK(tubeCylinderVolume->volumes().size() == 0u);
  BOOST_CHECK(tubeCylinderVolume->portals().size() == 4u);
}

BOOST_AUTO_TEST_SUITE_END()
