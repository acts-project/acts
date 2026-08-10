// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <stdexcept>
#include <variant>

using Acts::AxisFactory;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(AxisFactoryEquidistant) {
  using enum Acts::AxisBoundaryType;

  AxisFactory af = AxisFactory::Equidistant(5, 0., 10., Bound);
  BOOST_CHECK(!af.isDeferred());
  BOOST_CHECK(af.isEquidistant());
  BOOST_CHECK(!af.isVariable());
  BOOST_CHECK(!af.direction().has_value());
  BOOST_CHECK(af.boundaryType() == Bound);
  BOOST_CHECK_EQUAL(af.nBins(), 5);
  BOOST_CHECK(af.asEquidistant().min == 0.);
  BOOST_CHECK(af.asEquidistant().max == 10.);
  BOOST_CHECK_THROW(af.asVariable(), std::bad_variant_access);

  auto axis = af.buildAxis({});
  BOOST_CHECK(axis->isEquidistant());
  BOOST_CHECK_EQUAL(axis->getBoundaryType(), Bound);
  BOOST_CHECK_EQUAL(axis->getNBins(), 5);
  CHECK_CLOSE_ABS(axis->getMin(), 0., 1e-15);
  CHECK_CLOSE_ABS(axis->getMax(), 10., 1e-15);
  BOOST_CHECK(!axis->getDirection().has_value());

  // Supplying a described property validates it instead of overriding it
  BOOST_CHECK_NO_THROW(af.buildAxis({.min = 0., .max = 10.}));
  BOOST_CHECK_THROW(af.buildAxis({.min = 0., .max = 1.}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(af.buildAxis({.boundaryType = Closed}),
                    std::invalid_argument);

  // Invalid construction
  BOOST_CHECK_THROW(AxisFactory::Equidistant(5, 1., 0., Bound),
                    std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::Equidistant(0, 0., 1., Bound),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AxisFactoryVariable) {
  using enum Acts::AxisBoundaryType;

  AxisFactory af = AxisFactory::Variable({0., 1., 4., 10.}, Open);
  BOOST_CHECK(!af.isDeferred());
  BOOST_CHECK(!af.isEquidistant());
  BOOST_CHECK(af.isVariable());
  BOOST_CHECK(af.boundaryType() == Open);
  BOOST_CHECK_EQUAL(af.nBins(), 3);

  auto axis = af.buildAxis({});
  BOOST_CHECK(axis->isVariable());
  BOOST_CHECK_EQUAL(axis->getBoundaryType(), Open);
  std::vector<double> expectedEdges = {0., 1., 4., 10.};
  BOOST_CHECK(axis->getBinEdges() == expectedEdges);

  // The edges are the range, so supplying a different one is a mismatch
  BOOST_CHECK_THROW(af.buildAxis({.min = 0., .max = 1.}),
                    std::invalid_argument);

  // A variable axis without a boundary type takes it from the consumer
  AxisFactory afOpen = AxisFactory::Variable({0., 1., 4., 10.});
  BOOST_CHECK(afOpen.isDeferred());
  BOOST_CHECK_THROW(afOpen.buildAxis({}), std::domain_error);
  BOOST_CHECK_EQUAL(
      afOpen.buildAxis({.boundaryType = Bound})->getBoundaryType(), Bound);

  // Invalid construction
  BOOST_CHECK_THROW(AxisFactory::Variable({0.}, Bound), std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::Variable({0., 1., 1.}, Bound),
                    std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::Variable({0., 2., 1.}, Bound),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AxisFactoryDeferredEquidistant) {
  using enum Acts::AxisBoundaryType;

  AxisFactory af = AxisFactory::DeferredEquidistant(20);
  BOOST_CHECK(af.isDeferred());
  BOOST_CHECK(af.isEquidistant());
  BOOST_CHECK(!af.boundaryType().has_value());
  BOOST_CHECK_EQUAL(af.nBins(), 20);
  BOOST_CHECK_EQUAL(af.asEquidistant().nBins, 20);

  // A deferred description cannot be built without the missing properties
  BOOST_CHECK_THROW(af.buildAxis({}), std::domain_error);
  BOOST_CHECK_THROW(af.buildAxis({.min = -2., .max = 2.}), std::domain_error);

  auto axis = af.buildAxis({.min = -2., .max = 2., .boundaryType = Closed});
  BOOST_CHECK(axis->isEquidistant());
  BOOST_CHECK_EQUAL(axis->getBoundaryType(), Closed);
  BOOST_CHECK_EQUAL(axis->getNBins(), 20);
  CHECK_CLOSE_ABS(axis->getMin(), -2., 1e-15);
  CHECK_CLOSE_ABS(axis->getMax(), 2., 1e-15);

  // Invalid supplied range
  BOOST_CHECK_THROW(
      af.buildAxis({.min = 2., .max = -2., .boundaryType = Bound}),
      std::invalid_argument);

  // A range fixed at configuration time, the boundary type left open
  AxisFactory afRanged = AxisFactory::Equidistant(4, 0., 8.);
  BOOST_CHECK(afRanged.isDeferred());
  auto ranged = afRanged.buildAxis({.boundaryType = Bound});
  CHECK_CLOSE_ABS(ranged->getMin(), 0., 1e-15);
  CHECK_CLOSE_ABS(ranged->getMax(), 8., 1e-15);

  // Invalid construction
  BOOST_CHECK_THROW(AxisFactory::DeferredEquidistant(0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AxisFactoryDeferredVariable) {
  using enum Acts::AxisBoundaryType;

  AxisFactory af = AxisFactory::DeferredVariable({0., 0.1, 0.5, 1.});
  BOOST_CHECK(af.isDeferred());
  BOOST_CHECK(af.isVariable());
  BOOST_CHECK_EQUAL(af.nBins(), 3);

  BOOST_CHECK(af.isDeferredVariable());
  BOOST_CHECK_THROW(af.buildAxis({}), std::domain_error);

  // Edges are scaled affinely onto the supplied range
  auto axis = af.buildAxis({.min = 10., .max = 30., .boundaryType = Bound});
  BOOST_CHECK(axis->isVariable());
  auto edges = axis->getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 4);
  CHECK_CLOSE_ABS(edges[0], 10., 1e-15);
  CHECK_CLOSE_ABS(edges[1], 12., 1e-15);
  CHECK_CLOSE_ABS(edges[2], 20., 1e-15);
  CHECK_CLOSE_ABS(edges[3], 30., 1e-15);

  // Invalid construction: not normalized to [0, 1] or not strictly increasing
  BOOST_CHECK_THROW(AxisFactory::DeferredVariable({0., 0.5}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::DeferredVariable({0.1, 1.}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::DeferredVariable({0., 0.5, 0.4, 1.}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(AxisFactory::DeferredVariable({1.}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AxisFactoryFromAxis) {
  using enum Acts::AxisBoundaryType;

  auto eqAxis = Acts::IAxis::createEquidistant(Bound, -1., 1., 8,
                                               Acts::AxisDirection::AxisZ);
  AxisFactory af = AxisFactory::FromAxis(*eqAxis);
  BOOST_CHECK(!af.isDeferred());
  BOOST_CHECK(af.isEquidistant());
  BOOST_CHECK(af.direction() == Acts::AxisDirection::AxisZ);
  // Round trip
  BOOST_CHECK(*af.buildAxis({}) == *eqAxis);

  auto varAxis = Acts::IAxis::createVariable(Closed, {0., 2., 3.});
  AxisFactory afVar = AxisFactory::FromAxis(*varAxis);
  BOOST_CHECK(afVar.isVariable());
  BOOST_CHECK(!afVar.direction().has_value());
  BOOST_CHECK(*afVar.buildAxis({}) == *varAxis);
}

BOOST_AUTO_TEST_CASE(AxisFactoryDirectionHandling) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;

  AxisFactory af = AxisFactory::DeferredEquidistant(10, AxisPhi);
  BOOST_CHECK(af.direction() == AxisPhi);

  // Matching supplied direction is fine
  auto axis = af.buildAxis(
      {.min = 0., .max = 1., .boundaryType = Bound, .direction = AxisPhi});
  BOOST_CHECK(axis->getDirection() == AxisPhi);

  // Mismatching supplied direction throws
  BOOST_CHECK_THROW(
      af.buildAxis(
          {.min = 0., .max = 1., .boundaryType = Bound, .direction = AxisZ}),
      std::invalid_argument);

  // Without a described direction the supplied one is adopted
  AxisFactory afFree = AxisFactory::DeferredEquidistant(10);
  auto axisFree = afFree.buildAxis(
      {.min = 0., .max = 1., .boundaryType = Bound, .direction = AxisZ});
  BOOST_CHECK(axisFree->getDirection() == AxisZ);

  // withDirection attaches the direction
  AxisFactory afDir = afFree.withDirection(AxisRPhi);
  BOOST_CHECK(afDir.direction() == AxisRPhi);
  BOOST_CHECK(afFree != afDir);

  // Same rule for a description that leaves nothing else open
  AxisFactory afFull = AxisFactory::Equidistant(4, 0., 1., Bound, AxisX);
  BOOST_CHECK_THROW(afFull.buildAxis({.direction = AxisY}),
                    std::invalid_argument);
  BOOST_CHECK(afFull.buildAxis({.direction = AxisX})->getDirection() == AxisX);
  BOOST_CHECK(afFull.buildAxis({})->getDirection() == AxisX);
}

BOOST_AUTO_TEST_CASE(AxisFactoryToDeferred) {
  using enum Acts::AxisBoundaryType;

  // Equidistant keeps only the bin count
  AxisFactory af = AxisFactory::Equidistant(12, -3., 3., Closed,
                                            Acts::AxisDirection::AxisPhi);
  AxisFactory deferred = af.toDeferred();
  BOOST_CHECK(deferred.isDeferred());
  BOOST_CHECK(deferred.isEquidistant());
  BOOST_CHECK_EQUAL(deferred.nBins(), 12);
  BOOST_CHECK(deferred.direction() == Acts::AxisDirection::AxisPhi);

  // Variable edges are normalized to [0, 1] with exact endpoints
  AxisFactory afVar = AxisFactory::Variable({10., 12., 20., 30.}, Bound);
  AxisFactory deferredVar = afVar.toDeferred();
  BOOST_CHECK(deferredVar.isDeferred());
  const auto& normalizedEdges =
      deferredVar.asDeferredVariable().normalizedEdges;
  BOOST_CHECK_EQUAL(normalizedEdges.size(), 4);
  BOOST_CHECK_EQUAL(normalizedEdges.front(), 0.);
  BOOST_CHECK_EQUAL(normalizedEdges.back(), 1.);
  CHECK_CLOSE_ABS(normalizedEdges[1], 0.1, 1e-15);
  CHECK_CLOSE_ABS(normalizedEdges[2], 0.5, 1e-15);

  // Deferred descriptions are returned unchanged
  AxisFactory afDef = AxisFactory::DeferredEquidistant(7);
  BOOST_CHECK(afDef.toDeferred() == afDef);
}

BOOST_AUTO_TEST_CASE(AxisFactoryEqualityAndStreams) {
  using enum Acts::AxisBoundaryType;

  AxisFactory a = AxisFactory::Equidistant(10, 0., 1., Bound);
  AxisFactory b = AxisFactory::Equidistant(10, 0., 1., Bound);
  AxisFactory c = AxisFactory::Equidistant(10, 0., 2., Bound);
  BOOST_CHECK(a == b);
  BOOST_CHECK(a != c);
  BOOST_CHECK(a != AxisFactory::DeferredEquidistant(10));
  BOOST_CHECK(a != a.withDirection(Acts::AxisDirection::AxisX));

  BOOST_CHECK_EQUAL(a.toString(),
                    "AxisFactory: 10 bins, equidistant within [0, 1], Bound");
  BOOST_CHECK_EQUAL(
      AxisFactory::DeferredEquidistant(5, Acts::AxisDirection::AxisZ)
          .toString(),
      "AxisFactory: 5 bins in AxisZ, equidistant within deferred range, "
      "deferred boundary type");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
