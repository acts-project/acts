// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

using namespace Acts;
using namespace Acts::detail::Test;
using SourceLink = Acts::detail::Test::TestSourceLink;
namespace bd = boost::unit_test::data;

namespace {
constexpr BoundIndices boundIndices[] = {
    eBoundLoc0, eBoundLoc1, eBoundTime, eBoundPhi, eBoundTheta, eBoundQOverP,
};
constexpr FreeIndices freeIndices[] = {
    eFreePos0, eFreePos1, eFreePos2, eFreeTime,
    eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP,
};
const TestSourceLink sourceOrig;
const Acts::SourceLink source{sourceOrig};
// fix seed for reproducible tests
std::default_random_engine rng(123);
}  // namespace

// the underlying subspace implementation is already tested exhaustively in a
// separate unit test. here we only test concrete extreme cases and
// measurement-specific functionality.

BOOST_AUTO_TEST_SUITE(EventDataMeasurement)

BOOST_DATA_TEST_CASE(FixedBoundOne, bd::make(boundIndices), index) {
  auto [params, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
  auto meas = makeMeasurement(source, params, cov, index);

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : boundIndices) {
    if (i == index) {
      BOOST_CHECK(meas.contains(i));
    } else {
      BOOST_CHECK(!meas.contains(i));
    }
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink().template get<TestSourceLink>(),
                    sourceOrig);
}

BOOST_AUTO_TEST_CASE(FixedBoundAll) {
  auto [params, cov] = generateBoundParametersCovariance(rng);
  auto meas = makeMeasurement(source, params, cov, eBoundLoc0, eBoundLoc1,
                              eBoundPhi, eBoundTheta, eBoundQOverP, eBoundTime);

  BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
  for (auto i : boundIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink().get<TestSourceLink>(), sourceOrig);
}

namespace {
// example data for phi residual tests. each entry contains
//
//     measured, reference, expected residual
//
const std::vector<std::tuple<double, double, double>> kPhiDataset = {
    // measurement and reference in bounds and close
    {0.5, 0.75, -0.25},
    // measurement and reference in bounds but at different edges
    {0.25, 2 * M_PI - 0.25, 0.5},
    {2 * M_PI - 0.125, 0.125, -0.25},
    // measurement in bounds, reference ouf-of-bounds, both near lower edge
    {0.25, -0.25, 0.5},
    // measurement in bounds, reference ouf-of-bounds, both near upper edge
    {2 * M_PI - 0.25, 2 * M_PI + 0.25, -0.5},
    // measurement out-of-bounds, reference in bounds, both near lower edge
    {-0.25, 0.25, -0.5},
    // measurement out-of-bounds, reference in bounds, both near upper edge
    {2 * M_PI + 0.25, 2 * M_PI - 0.25, 0.5},
};
}  // namespace

BOOST_DATA_TEST_CASE(BoundResidualsPhi, bd::make(kPhiDataset), phiMea, phiRef,
                     phiRes) {
  using MeasurementVector = Acts::ActsVector<1>;
  using MeasurementCovariance = Acts::ActsSquareMatrix<1>;

  // prepare measurement
  MeasurementVector params = MeasurementVector::Zero();
  MeasurementCovariance cov = MeasurementCovariance::Zero();
  params[0] = phiMea;
  auto measurement = makeMeasurement(source, params, cov, eBoundPhi);
  // prepare reference parameters
  Acts::BoundVector reference = Acts::BoundVector::Zero();
  reference[eBoundPhi] = phiRef;

  // compute and check residual
  auto res = measurement.residuals(reference);
  CHECK_CLOSE_ABS(res[0], phiRes, std::numeric_limits<ActsScalar>::epsilon());
}

BOOST_DATA_TEST_CASE(FixedFreeOne, bd::make(freeIndices), index) {
  auto [params, cov] = generateParametersCovariance<ActsScalar, 1u>(rng);
  auto meas = makeMeasurement(source, params, cov, index);

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : freeIndices) {
    if (i == index) {
      BOOST_CHECK(meas.contains(i));
    } else {
      BOOST_CHECK(!meas.contains(i));
    }
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink().template get<TestSourceLink>(),
                    sourceOrig);

  // all free parameters are unrestricted and we know the expected residual.
  constexpr auto tol = std::numeric_limits<ActsScalar>::epsilon();
  auto [ref, refCov] = generateFreeParametersCovariance(rng);
  auto res = meas.residuals(ref);
  CHECK_CLOSE_ABS(res[0], params[0] - ref[index], tol);
}

BOOST_AUTO_TEST_CASE(FixedFreeAll) {
  auto [params, cov] = generateFreeParametersCovariance(rng);
  auto meas =
      makeMeasurement(source, params, cov, eFreePos0, eFreePos1, eFreePos2,
                      eFreeTime, eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP);

  BOOST_CHECK_EQUAL(meas.size(), eFreeSize);
  for (auto i : freeIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink().get<TestSourceLink>(), sourceOrig);

  // all free parameters are unrestricted and we know the expected residual.
  constexpr auto tol = std::numeric_limits<ActsScalar>::epsilon();
  auto [ref, refCov] = generateFreeParametersCovariance(rng);
  CHECK_CLOSE_ABS(meas.residuals(ref), params - ref, tol);
}

BOOST_AUTO_TEST_CASE(VariantBound) {
  // generate w/ a single parameter
  auto [par1, cov1] = generateParametersCovariance<ActsScalar, 1u>(rng);
  BoundVariantMeasurement meas =
      makeMeasurement(source, par1, cov1, eBoundTheta);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), 1);
        BOOST_CHECK(!m.contains(eBoundLoc0));
        BOOST_CHECK(!m.contains(eBoundLoc1));
        BOOST_CHECK(!m.contains(eBoundTime));
        BOOST_CHECK(!m.contains(eBoundPhi));
        BOOST_CHECK(m.contains(eBoundTheta));
        BOOST_CHECK(!m.contains(eBoundQOverP));
      },
      meas);

  // reassign w/ all parameters
  auto [parN, covN] = generateBoundParametersCovariance(rng);
  meas = makeMeasurement(source, parN, covN, eBoundLoc0, eBoundLoc1, eBoundPhi,
                         eBoundTheta, eBoundQOverP, eBoundTime);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), eBoundSize);
        BOOST_CHECK(m.contains(eBoundLoc0));
        BOOST_CHECK(m.contains(eBoundLoc1));
        BOOST_CHECK(m.contains(eBoundTime));
        BOOST_CHECK(m.contains(eBoundPhi));
        BOOST_CHECK(m.contains(eBoundTheta));
        BOOST_CHECK(m.contains(eBoundQOverP));
      },
      meas);
}

BOOST_AUTO_TEST_CASE(VariantFree) {
  // generate w/ two parameters
  auto [par2, cov2] = generateParametersCovariance<ActsScalar, 2u>(rng);
  FreeVariantMeasurement meas =
      makeMeasurement(source, par2, cov2, eFreePos2, eFreeTime);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), 2);
        BOOST_CHECK(!m.contains(eFreePos0));
        BOOST_CHECK(!m.contains(eFreePos1));
        BOOST_CHECK(m.contains(eFreePos2));
        BOOST_CHECK(m.contains(eFreeTime));
        BOOST_CHECK(!m.contains(eFreeDir0));
        BOOST_CHECK(!m.contains(eFreeDir1));
        BOOST_CHECK(!m.contains(eFreeDir2));
        BOOST_CHECK(!m.contains(eFreeQOverP));
      },
      meas);

  // reassign w/ all parameters
  auto [parN, covN] = generateFreeParametersCovariance(rng);
  meas =
      makeMeasurement(source, parN, covN, eFreePos0, eFreePos1, eFreePos2,
                      eFreeTime, eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), eFreeSize);
        BOOST_CHECK(m.contains(eFreePos0));
        BOOST_CHECK(m.contains(eFreePos1));
        BOOST_CHECK(m.contains(eFreePos2));
        BOOST_CHECK(m.contains(eFreeTime));
        BOOST_CHECK(m.contains(eFreeDir0));
        BOOST_CHECK(m.contains(eFreeDir1));
        BOOST_CHECK(m.contains(eFreeDir2));
        BOOST_CHECK(m.contains(eFreeQOverP));
      },
      meas);
}

BOOST_AUTO_TEST_SUITE_END()
