// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"

#include <random>

namespace {

using namespace Acts;
using namespace Acts::Test;
using SourceLink = Acts::Test::TestSourceLink;
namespace bd = boost::unit_test::data;

constexpr BoundIndices boundIndices[] = {
    eBoundLoc0, eBoundLoc1, eBoundTime, eBoundPhi, eBoundTheta, eBoundQOverP,
};
constexpr FreeIndices freeIndices[] = {
    eFreePos0, eFreePos1, eFreePos2, eFreeTime,
    eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP,
};
const TestSourceLink source;
// fix seed for reproducible tests
std::default_random_engine rng(123);

}  // namespace

// the underlying subspace implementation is already tested exhaustively in a
// separate unit test. here we only test concrete extreme cases and
// measurement-specific functionality.

BOOST_AUTO_TEST_SUITE(EventDataMeasurement)

BOOST_DATA_TEST_CASE(FixedBoundOne, bd::make(boundIndices), index) {
  auto [params, cov] = generateParametersCovariance(rng, index);
  auto meas = makeFixedSizeMeasurement(source, params, cov, index);

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : boundIndices) {
    if (i == index) {
      BOOST_CHECK(meas.contains(i));
    } else {
      BOOST_CHECK(not meas.contains(i));
    }
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink(), source);
}

BOOST_AUTO_TEST_CASE(FixedBoundAll) {
  auto [params, cov] = generateBoundParametersCovariance(rng);
  auto meas = makeFixedSizeMeasurement(source, params, cov, eBoundLoc0,
                                       eBoundLoc1, eBoundPhi, eBoundTheta,
                                       eBoundQOverP, eBoundTime);

  BOOST_CHECK_EQUAL(meas.size(), eBoundSize);
  for (auto i : boundIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink(), source);
}

BOOST_DATA_TEST_CASE(FixedFreeOne, bd::make(freeIndices), index) {
  auto [params, cov] = generateParametersCovariance(rng, index);
  auto meas = makeFixedSizeMeasurement(source, params, cov, index);

  BOOST_CHECK_EQUAL(meas.size(), 1);
  for (auto i : freeIndices) {
    if (i == index) {
      BOOST_CHECK(meas.contains(i));
    } else {
      BOOST_CHECK(not meas.contains(i));
    }
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink(), source);
}

BOOST_AUTO_TEST_CASE(FixedFreeAll) {
  auto [params, cov] = generateFreeParametersCovariance(rng);
  auto meas = makeFixedSizeMeasurement(
      source, params, cov, eFreePos0, eFreePos1, eFreePos2, eFreeTime,
      eFreeDir0, eFreeDir1, eFreeDir2, eFreeQOverP);

  BOOST_CHECK_EQUAL(meas.size(), eFreeSize);
  for (auto i : freeIndices) {
    BOOST_CHECK(meas.contains(i));
  }
  BOOST_CHECK_EQUAL(meas.parameters(), params);
  BOOST_CHECK_EQUAL(meas.covariance(), cov);
  BOOST_CHECK_EQUAL(meas.sourceLink(), source);
}

BOOST_AUTO_TEST_CASE(VariantBound) {
  // generate w/ a single parameter
  auto [par1, cov1] = generateParametersCovariance(rng, eBoundTheta);
  auto meas = makeVariantMeasurement(source, par1, cov1, eBoundTheta);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), 1);
        BOOST_CHECK(not m.contains(eBoundLoc0));
        BOOST_CHECK(not m.contains(eBoundLoc1));
        BOOST_CHECK(not m.contains(eBoundTime));
        BOOST_CHECK(not m.contains(eBoundPhi));
        BOOST_CHECK(m.contains(eBoundTheta));
        BOOST_CHECK(not m.contains(eBoundQOverP));
      },
      meas);

  // reassign w/ all parameters
  auto [parN, covN] = generateBoundParametersCovariance(rng);
  meas =
      makeVariantMeasurement(source, parN, covN, eBoundLoc0, eBoundLoc1,
                             eBoundPhi, eBoundTheta, eBoundQOverP, eBoundTime);
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
  auto [par2, cov2] = generateParametersCovariance(rng, eFreePos2, eFreeTime);
  auto meas = makeVariantMeasurement(source, par2, cov2, eFreePos2, eFreeTime);
  std::visit(
      [](const auto& m) {
        BOOST_CHECK_EQUAL(m.size(), 2);
        BOOST_CHECK(not m.contains(eFreePos0));
        BOOST_CHECK(not m.contains(eFreePos1));
        BOOST_CHECK(m.contains(eFreePos2));
        BOOST_CHECK(m.contains(eFreeTime));
        BOOST_CHECK(not m.contains(eFreeDir0));
        BOOST_CHECK(not m.contains(eFreeDir1));
        BOOST_CHECK(not m.contains(eFreeDir2));
        BOOST_CHECK(not m.contains(eFreeQOverP));
      },
      meas);

  // reassign w/ all parameters
  auto [parN, covN] = generateFreeParametersCovariance(rng);
  meas = makeVariantMeasurement(source, parN, covN, eFreePos0, eFreePos1,
                                eFreePos2, eFreeTime, eFreeDir0, eFreeDir1,
                                eFreeDir2, eFreeQOverP);
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
