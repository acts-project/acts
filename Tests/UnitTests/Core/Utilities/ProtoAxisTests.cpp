// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

BOOST_AUTO_TEST_SUITE(ProtoAxis)

BOOST_AUTO_TEST_CASE(EquidistantProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis epab(AxisX, Bound, 0.0, 1.0, 10);

  BOOST_CHECK(epab.isEquidistant());

  BOOST_CHECK(!epab.isVariable());

  auto edges = epab.getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 11);

  BOOST_CHECK_EQUAL(epab.getAxisDirection(), AxisX);

  BOOST_CHECK_EQUAL(epab.getType(), Equidistant);

  BOOST_CHECK_EQUAL(epab.getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epab.getNBins(), 10);

  CHECK_CLOSE_ABS(epab.getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(epab.getMax(), 1.0, 1e-15);

  BOOST_CHECK(!epab.isAutorange());

  // Generate a 1 dimensional grid from double,
  // return type needs to be type-erased

  // Open, equidistant axis
  Acts::ProtoAxis epao(AxisY, Open, 0., 2.0, 10.);

  BOOST_CHECK_EQUAL(epao.getBoundaryType(), Open);

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, 1.0, 0.0, 10),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, 0.0, 1.0, 0),
                    std::invalid_argument);

  // Invalid constructor, closed with somethgin else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, 0.0, 1.0, 10),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AutorangeProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis, autorange
  Acts::ProtoAxis epa(AxisX, Bound, 10);

  BOOST_CHECK(epa.isEquidistant());

  BOOST_CHECK(!epa.isVariable());

  BOOST_CHECK_EQUAL(epa.getAxisDirection(), AxisX);

  BOOST_CHECK_EQUAL(epa.getType(), Equidistant);

  BOOST_CHECK_EQUAL(epa.getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epa.getNBins(), 10);

  BOOST_CHECK(epa.isAutorange());

  // Invalid constructor, closed with somethgin else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, 10), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VariabletProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis vpab(AxisX, Bound, {0.0, 1.0, 10});

  BOOST_CHECK(!vpab.isEquidistant());

  BOOST_CHECK(vpab.isVariable());

  auto edges = vpab.getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 3);

  BOOST_CHECK_EQUAL(vpab.getAxisDirection(), AxisX);

  BOOST_CHECK_EQUAL(vpab.getType(), Variable);

  BOOST_CHECK_EQUAL(vpab.getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(vpab.getNBins(), 2);

  CHECK_CLOSE_ABS(vpab.getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(vpab.getMax(), 10.0, 1e-15);

  BOOST_CHECK(!vpab.isAutorange());

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, std::vector<double>{2.}),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Bound, {3., 2., 1}),
                    std::invalid_argument);

  // Invalid constructor, closed with somethgin else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, {0., 1., 2., 3.}),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
