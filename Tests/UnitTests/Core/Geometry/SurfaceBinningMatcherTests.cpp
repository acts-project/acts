// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/SurfaceBinningMatcher.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <cmath>
#include <memory>
#include <vector>

#include <boost/format.hpp>

namespace Acts::Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_CASE(PlaneSurfaceMatcher) {
  auto identity = Transform3::Identity();

  double rMin = 5.;
  double rMax = 10.;
  double rMinTol = 0.1;
  double rMaxTol = 0.5;

  double phiTol = 0.1;

  auto oneBounds = std::make_shared<RadialBounds>(rMin, rMax, M_PI / 16, 0.);
  auto oneSurface = Surface::makeShared<DiscSurface>(identity, oneBounds);

  auto otherBounds =
      std::make_shared<RadialBounds>(2 * rMax, 4 * rMax, M_PI / 16, 0.5 * M_PI);
  auto otherSurface = Surface::makeShared<DiscSurface>(identity, otherBounds);

  auto similarRbounds = std::make_shared<RadialBounds>(
      rMin - 0.5 * rMinTol, rMax + 0.5 * rMaxTol, M_PI / 16, 0.5 * M_PI);
  auto similarRSurface =
      Surface::makeShared<DiscSurface>(identity, similarRbounds);

  auto similarPhiBounds =
      std::make_shared<RadialBounds>(0.25 * rMin, 0.5 * rMin, M_PI / 16, 0.);
  auto similarPhiSurface =
      Surface::makeShared<DiscSurface>(identity, similarPhiBounds);

  SurfaceBinningMatcher sbm;
  sbm.tolerances[binR] = {rMinTol, rMaxTol};
  sbm.tolerances[binPhi] = {phiTol, phiTol};

  // Always true
  for (int ib = 0; ib < binValues; ++ib) {
    BOOST_CHECK(
        sbm(tgContext, (BinningValue)ib, oneSurface.get(), oneSurface.get()));
  }
  // Not matching in R
  BOOST_CHECK(!sbm(tgContext, binR, oneSurface.get(), otherSurface.get()));
  // Not matching in phi
  BOOST_CHECK(!sbm(tgContext, binPhi, oneSurface.get(), otherSurface.get()));

  // Good enough matching in R
  BOOST_CHECK(sbm(tgContext, binR, oneSurface.get(), similarRSurface.get()));
  // Good enough matching in phi
  BOOST_CHECK(
      sbm(tgContext, binPhi, oneSurface.get(), similarPhiSurface.get()));
}

}  // namespace Acts::Test
