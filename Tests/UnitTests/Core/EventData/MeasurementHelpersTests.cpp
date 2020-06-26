// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

using SourceLink = MinimalSourceLink;

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;
using FittableMeasurement = FittableMeasurement<SourceLink>;

BOOST_AUTO_TEST_CASE(getSurface_test) {
  auto cylinderBounds = std::make_shared<CylinderBounds>(3, 10);

  auto cylinder = Surface::makeShared<CylinderSurface>(nullptr, cylinderBounds);
  auto cylinder2 =
      Surface::makeShared<CylinderSurface>(nullptr, cylinderBounds);

  SymMatrix2D cov;
  cov << 0.04, 0, 0, 0.1;
  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m(cylinder, {},
                                                    std::move(cov), -0.1, 0.45);

  FittableMeasurement fm = m;

  BOOST_CHECK_EQUAL(MeasurementHelpers::getSurface(fm), cylinder.get());

  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m2(
      cylinder2, {}, std::move(cov), -0.1, 0.45);
  fm = m2;
  BOOST_CHECK_EQUAL(MeasurementHelpers::getSurface(fm), cylinder2.get());
}

BOOST_AUTO_TEST_CASE(getSize_test) {
  auto cylinder = Surface::makeShared<CylinderSurface>(nullptr, 3, 10);

  SymMatrix2D cov;
  cov << 0.04, 0, 0, 0.1;
  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m(cylinder, {},
                                                    std::move(cov), -0.1, 0.45);

  FittableMeasurement fm = m;
  BOOST_CHECK_EQUAL(MeasurementHelpers::getSize(fm), 2u);

  ActsSymMatrixD<3> cov3;
  cov.setRandom();
  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1, ParDef::eT> m2(
      cylinder, {}, std::move(cov3), -0.1, 0.45, 42);
  fm = m2;

  BOOST_CHECK_EQUAL(MeasurementHelpers::getSize(fm), 3u);
}

BOOST_AUTO_TEST_CASE(MinimalSourceLinkTest) {
  auto cylinder = Surface::makeShared<CylinderSurface>(nullptr, 3, 10);

  SymMatrix2D cov;
  cov << 0.04, 0, 0, 0.1;
  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m(cylinder, {},
                                                    std::move(cov), -0.1, 0.45);

  FittableMeasurement fm = m;
  MinimalSourceLink msl{&fm};

  BOOST_CHECK_EQUAL(&msl.referenceSurface(), cylinder.get());

  MinimalSourceLink msl2{&fm};
  BOOST_CHECK_EQUAL(msl, msl2);

  MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m2(
      cylinder, {}, std::move(cov), -0.1, 0.45);
  FittableMeasurement fm2 = m2;
  MinimalSourceLink msl3{&fm2};
  BOOST_CHECK_NE(msl, msl3);

  BOOST_CHECK_EQUAL(&*msl, &fm);
}

BOOST_AUTO_TEST_CASE(visit_measurement_test) {
  // Overallocated full size parameter vector and covariance
  BoundVector parFull = BoundVector::Random();
  BoundMatrix covFull = BoundMatrix::Random();
  // constant variants
  const auto& parFullConst = parFull;
  const auto& covFullConst = covFull;

  for (BoundVector::Index dim = 1; dim <= parFull.size(); ++dim) {
    visit_measurement(parFull, covFull, dim, [&](auto param, auto cov) {
      BOOST_CHECK_EQUAL(param, parFull.head(dim));
      BOOST_CHECK_EQUAL(cov, covFull.topLeftCorner(dim, dim));
    });
    visit_measurement(parFull, covFull, dim,
                      [&](const auto& param, const auto& cov) {
                        BOOST_CHECK_EQUAL(param, parFull.head(dim));
                        BOOST_CHECK_EQUAL(cov, covFull.topLeftCorner(dim, dim));
                      });
    visit_measurement(parFullConst, covFullConst, dim,
                      [&](const auto& param, const auto& cov) {
                        BOOST_CHECK_EQUAL(param, parFullConst.head(dim));
                        BOOST_CHECK_EQUAL(cov,
                                          covFullConst.topLeftCorner(dim, dim));
                      });
  }
}

}  // namespace Test
}  // namespace Acts
