// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/MinimalSourceLink.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

BOOST_AUTO_TEST_SUITE(EventDataMinimalSourceLink)

BOOST_AUTO_TEST_CASE(ConstructAndAssign) {
  using namespace Acts;
  using ThisMeasurement =
      Measurement<MinimalSourceLink, BoundIndices, eBoundLoc0, eBoundLoc1>;
  using ThisFittableMeasurement = FittableMeasurement<MinimalSourceLink>;

  auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3D::Identity(), 3, 10);

  SymMatrix2D cov;
  cov << 0.04, 0, 0, 0.1;

  ThisMeasurement m(cylinder, {}, std::move(cov), -0.1, 0.45);

  ThisFittableMeasurement fm = m;
  MinimalSourceLink msl{&fm};

  BOOST_CHECK_EQUAL(msl.geometryId(), cylinder->geometryId());

  MinimalSourceLink msl2{&fm};
  BOOST_CHECK_EQUAL(msl, msl2);

  ThisMeasurement m2(cylinder, {}, std::move(cov), -0.1, 0.45);
  ThisFittableMeasurement fm2 = m2;
  MinimalSourceLink msl3{&fm2};
  BOOST_CHECK_NE(msl, msl3);

  BOOST_CHECK_EQUAL(msl.measurement, &fm);
}

BOOST_AUTO_TEST_SUITE_END()
