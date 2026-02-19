// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include <ActsPodioEdm/TrackerHitLocalCollection.h>

using namespace Acts;
using namespace ActsPlugins;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::Test;

namespace {
std::default_random_engine rng(123);
auto gctx = GeometryContext::dangerouslyDefaultConstruct();
}  // namespace

BOOST_AUTO_TEST_SUITE(EDM4hepMeasurementTest)

BOOST_AUTO_TEST_CASE(WriteMeasurement) {
  auto [parameters, covariance] =
      generateParametersCovariance<double, eBoundSize>(rng);

  std::uint64_t cellId = 1234;

  ActsPodioEdm::TrackerHitLocalCollection hits;
  auto to = hits.create();

  std::vector<std::uint8_t> indices = {eBoundLoc0, eBoundLoc1};
  FixedSubspaceHelper<eBoundSize, 2> helper(indices);

  Vector2 measPos = helper.projectVector(parameters);
  SquareMatrix<2> measCov = helper.projectMatrix(covariance);

  auto surface = Surface::makeShared<DiscSurface>(
      Transform3::Identity() * Translation3(Vector3{1, 2, 3}), 1_mm, 1_mm);

  Vector3 global = surface->localToGlobal(gctx, measPos, Vector3::UnitZ());

  EDM4hepUtil::writeMeasurement(gctx, {measPos.data(), 2},
                                {measCov.data(), 2, 2}, indices, cellId,
                                *surface, to);

  BOOST_CHECK_EQUAL(to.getCellID(), cellId);
  BOOST_CHECK_EQUAL(to.getPosition()[0] * 1_mm, global.x());
  BOOST_CHECK_EQUAL(to.getPosition()[1] * 1_mm, global.y());
  BOOST_CHECK_EQUAL(to.getPosition()[2] * 1_mm, global.z());

  // Time should be zero since we don't have time measurement
  BOOST_CHECK_EQUAL(to.getTime() * 1_ns, 0.0);

  auto meas = to.getMeasurement();
  BOOST_CHECK_EQUAL(meas.size(), 2);
  CHECK_CLOSE_REL(meas[0], measPos.x(), 1e-6);
  CHECK_CLOSE_REL(meas[1], measPos.y(), 1e-6);

  auto cov = to.getCovariance();
  BOOST_CHECK_EQUAL(cov.size(), 4);
  CHECK_CLOSE_REL(cov[0], measCov(0, 0), 1e-6);
  CHECK_CLOSE_REL(cov[1], measCov(0, 1), 1e-6);
  CHECK_CLOSE_REL(cov[2], measCov(1, 0), 1e-6);
  CHECK_CLOSE_REL(cov[3], measCov(1, 1), 1e-6);

  auto unpackedIndices = EDM4hepUtil::detail::decodeIndices(to.getType());
  BOOST_CHECK_EQUAL(unpackedIndices.size(), 2);
  BOOST_CHECK_EQUAL(unpackedIndices[0], eBoundLoc0);
  BOOST_CHECK_EQUAL(unpackedIndices[1], eBoundLoc1);

  // Round-trip: read back and verify
  auto read = EDM4hepUtil::readMeasurement(to);
  BOOST_CHECK_EQUAL(read.cellId, cellId);
  BOOST_CHECK_EQUAL(read.indices.size(), 2);
  BOOST_CHECK_EQUAL(read.indices[0], eBoundLoc0);
  BOOST_CHECK_EQUAL(read.indices[1], eBoundLoc1);
  BOOST_CHECK_EQUAL(read.parameters.size(), 2);
  CHECK_CLOSE_REL(read.parameters(0), measPos.x(), 1e-6);
  CHECK_CLOSE_REL(read.parameters(1), measPos.y(), 1e-6);
  BOOST_CHECK_EQUAL(read.covariance.rows(), 2);
  BOOST_CHECK_EQUAL(read.covariance.cols(), 2);
  CHECK_CLOSE_REL(read.covariance(0, 0), measCov(0, 0), 1e-6);
  CHECK_CLOSE_REL(read.covariance(0, 1), measCov(0, 1), 1e-6);
  CHECK_CLOSE_REL(read.covariance(1, 0), measCov(1, 0), 1e-6);
  CHECK_CLOSE_REL(read.covariance(1, 1), measCov(1, 1), 1e-6);
}

BOOST_AUTO_TEST_CASE(WriteMeasurementNoPosition) {
  auto [parameters, covariance] =
      generateParametersCovariance<double, eBoundSize>(rng);

  std::uint64_t cellId = 1234;

  ActsPodioEdm::TrackerHitLocalCollection hits;
  auto to = hits.create();

  // Only measure phi and theta
  std::vector<std::uint8_t> indices = {eBoundPhi, eBoundTheta};
  FixedSubspaceHelper<eBoundSize, 2> helper(indices);

  Vector2 measPos = helper.projectVector(parameters);
  SquareMatrix<2> measCov = helper.projectMatrix(covariance);

  auto surface = Surface::makeShared<DiscSurface>(
      Transform3::Identity() * Translation3(Vector3{1, 2, 3}), 1_mm, 1_mm);

  EDM4hepUtil::writeMeasurement(gctx, {measPos.data(), 2},
                                {measCov.data(), 2, 2}, indices, cellId,
                                *surface, to);

  BOOST_CHECK_EQUAL(to.getCellID(), cellId);
  // Position should be zero since we don't have loc0/loc1
  BOOST_CHECK_EQUAL(to.getPosition()[0] * 1_mm, 0.0);
  BOOST_CHECK_EQUAL(to.getPosition()[1] * 1_mm, 0.0);
  BOOST_CHECK_EQUAL(to.getPosition()[2] * 1_mm, 0.0);

  auto meas = to.getMeasurement();
  BOOST_CHECK_EQUAL(meas.size(), 2);
  CHECK_CLOSE_REL(meas[0], measPos.x(), 1e-6);
  CHECK_CLOSE_REL(meas[1], measPos.y(), 1e-6);

  auto cov = to.getCovariance();
  BOOST_CHECK_EQUAL(cov.size(), 4);
  CHECK_CLOSE_REL(cov[0], measCov(0, 0), 1e-6);
  CHECK_CLOSE_REL(cov[1], measCov(0, 1), 1e-6);
  CHECK_CLOSE_REL(cov[2], measCov(1, 0), 1e-6);
  CHECK_CLOSE_REL(cov[3], measCov(1, 1), 1e-6);

  auto unpackedIndices = EDM4hepUtil::detail::decodeIndices(to.getType());
  BOOST_CHECK_EQUAL(unpackedIndices.size(), 2);
  BOOST_CHECK_EQUAL(unpackedIndices[0], eBoundPhi);
  BOOST_CHECK_EQUAL(unpackedIndices[1], eBoundTheta);
}

BOOST_AUTO_TEST_CASE(WriteMeasurementWithTime) {
  auto [parameters, covariance] =
      generateParametersCovariance<double, eBoundSize>(rng);

  std::uint64_t cellId = 1234;

  ActsPodioEdm::TrackerHitLocalCollection hits;
  auto to = hits.create();

  // Measure loc0, loc1, and time
  std::vector<std::uint8_t> indices = {eBoundLoc0, eBoundLoc1, eBoundTime};
  FixedSubspaceHelper<eBoundSize, 3> helper(indices);

  Vector3 measPos = helper.projectVector(parameters);
  SquareMatrix<3> measCov = helper.projectMatrix(covariance);

  auto surface = Surface::makeShared<DiscSurface>(
      Transform3::Identity() * Translation3(Vector3{1, 2, 3}), 1_mm, 1_mm);

  Vector3 global =
      surface->localToGlobal(gctx, measPos.head<2>(), Vector3::UnitZ());

  EDM4hepUtil::writeMeasurement(gctx, {measPos.data(), 3},
                                {measCov.data(), 3, 3}, indices, cellId,
                                *surface, to);

  BOOST_CHECK_EQUAL(to.getCellID(), cellId);
  BOOST_CHECK_EQUAL(to.getPosition()[0] * 1_mm, global.x());
  BOOST_CHECK_EQUAL(to.getPosition()[1] * 1_mm, global.y());
  BOOST_CHECK_EQUAL(to.getPosition()[2] * 1_mm, global.z());
  // Time should be set since we have time measurement
  CHECK_CLOSE_REL(to.getTime() * 1_ns, parameters[eBoundTime], 1e-6);

  auto meas = to.getMeasurement();
  BOOST_CHECK_EQUAL(meas.size(), 3);
  CHECK_CLOSE_REL(meas[0], measPos.x(), 1e-6);
  CHECK_CLOSE_REL(meas[1], measPos.y(), 1e-6);
  CHECK_CLOSE_REL(meas[2], measPos.z(), 1e-6);

  auto cov = to.getCovariance();
  BOOST_CHECK_EQUAL(cov.size(), 9);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      CHECK_CLOSE_REL(cov[i * 3 + j], measCov(i, j), 1e-6);
    }
  }

  auto unpackedIndices = EDM4hepUtil::detail::decodeIndices(to.getType());
  BOOST_CHECK_EQUAL(unpackedIndices.size(), 3);
  BOOST_CHECK_EQUAL(unpackedIndices[0], eBoundLoc0);
  BOOST_CHECK_EQUAL(unpackedIndices[1], eBoundLoc1);
  BOOST_CHECK_EQUAL(unpackedIndices[2], eBoundTime);
}

BOOST_AUTO_TEST_CASE(EncodeDecodeIndices) {
  // Test empty span
  {
    std::vector<std::uint8_t> indices = {};
    std::uint32_t encoded = EDM4hepUtil::detail::encodeIndices(indices);
    auto decoded = EDM4hepUtil::detail::decodeIndices(encoded);
    BOOST_CHECK_EQUAL(decoded.size(), 0);
  }

  // Test single value
  {
    std::vector<std::uint8_t> indices = {3};
    std::uint32_t encoded = EDM4hepUtil::detail::encodeIndices(indices);
    auto decoded = EDM4hepUtil::detail::decodeIndices(encoded);
    BOOST_CHECK_EQUAL(decoded.size(), 1);
    BOOST_CHECK_EQUAL(decoded.at(0), 3);
  }

  // Test maximum length (6)
  {
    std::vector<std::uint8_t> indices = {0, 1, 2, 3, 4, 5};
    auto encoded = EDM4hepUtil::detail::encodeIndices(indices);
    auto decoded = EDM4hepUtil::detail::decodeIndices(encoded);
    BOOST_CHECK_EQUAL(decoded.size(), 6);
    for (std::size_t i = 0; i < 6; ++i) {
      BOOST_CHECK_EQUAL(decoded.at(i), i);
    }
  }

  // Test maximum value (6)
  {
    std::vector<std::uint8_t> indices = {6, 6, 6};
    auto encoded = EDM4hepUtil::detail::encodeIndices(indices);
    auto decoded = EDM4hepUtil::detail::decodeIndices(encoded);
    BOOST_CHECK_EQUAL(decoded.size(), 3);
    for (std::size_t i = 0; i < 3; ++i) {
      BOOST_CHECK_EQUAL(decoded.at(i), 6);
    }
  }

  // Test mixed values
  {
    std::vector<std::uint8_t> indices = {2, 5, 1, 4};
    auto encoded = EDM4hepUtil::detail::encodeIndices(indices);
    auto decoded = EDM4hepUtil::detail::decodeIndices(encoded);
    BOOST_CHECK_EQUAL(decoded.size(), 4);
    BOOST_CHECK_EQUAL(decoded.at(0), 2);
    BOOST_CHECK_EQUAL(decoded.at(1), 5);
    BOOST_CHECK_EQUAL(decoded.at(2), 1);
    BOOST_CHECK_EQUAL(decoded.at(3), 4);
  }
}

BOOST_AUTO_TEST_CASE(EncodeDecodeIndicesErrors) {
  // Test exceeding maximum length (7 values)
  {
    std::vector<std::uint8_t> indices = {0, 1, 2, 3, 4, 5, 6};
    BOOST_CHECK_THROW(EDM4hepUtil::detail::encodeIndices(indices),
                      std::runtime_error);
  }

  // Test exceeding maximum value (7)
  {
    std::vector<std::uint8_t> indices = {7};
    BOOST_CHECK_THROW(EDM4hepUtil::detail::encodeIndices(indices),
                      std::runtime_error);
  }

  // Test mixed valid/invalid values
  {
    std::vector<std::uint8_t> indices = {2, 7, 1};
    BOOST_CHECK_THROW(EDM4hepUtil::detail::encodeIndices(indices),
                      std::runtime_error);
  }
}

BOOST_AUTO_TEST_SUITE_END()
