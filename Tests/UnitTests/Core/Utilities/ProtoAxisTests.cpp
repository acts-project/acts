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

  // Direct access
  BOOST_CHECK_EQUAL(epab.getAxisDirection(), AxisX);
  BOOST_CHECK(!epab.isAutorange());

  // Access via IAxis
  BOOST_CHECK(epab.getAxis().isEquidistant());

  BOOST_CHECK(!epab.getAxis().isVariable());

  auto edges = epab.getAxis().getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 11);

  BOOST_CHECK_EQUAL(epab.getAxis().getType(), Equidistant);

  BOOST_CHECK_EQUAL(epab.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epab.getAxis().getNBins(), 10);

  CHECK_CLOSE_ABS(epab.getAxis().getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(epab.getAxis().getMax(), 1.0, 1e-15);

  // Open, equidistant axis
  Acts::ProtoAxis epao(AxisY, Open, 0., 2.0, 10.);

  BOOST_CHECK_EQUAL(epao.getAxis().getBoundaryType(), Open);

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

  // Direct access
  BOOST_CHECK_EQUAL(epa.getAxisDirection(), AxisX);

  BOOST_CHECK(epa.isAutorange());

  // Access via IAxis
  BOOST_CHECK(epa.getAxis().isEquidistant());

  BOOST_CHECK(!epa.getAxis().isVariable());

  BOOST_CHECK_EQUAL(epa.getAxis().getType(), Equidistant);

  BOOST_CHECK_EQUAL(epa.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epa.getAxis().getNBins(), 10);

  // Invalid constructor, closed with somethgin else than phi or rphi
  BOOST_CHECK_THROW(Acts::ProtoAxis(AxisZ, Closed, 10), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VariabletProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis vpab(AxisX, Bound, {0.0, 1.0, 10});

  // Direct access
  BOOST_CHECK_EQUAL(vpab.getAxisDirection(), AxisX);

  BOOST_CHECK(!vpab.isAutorange());

  // Access via IAxis
  BOOST_CHECK(!vpab.getAxis().isEquidistant());

  BOOST_CHECK(vpab.getAxis().isVariable());

  auto edges = vpab.getAxis().getBinEdges();
  BOOST_CHECK_EQUAL(edges.size(), 3);

  BOOST_CHECK_EQUAL(vpab.getAxis().getType(), Variable);

  BOOST_CHECK_EQUAL(vpab.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(vpab.getAxis().getNBins(), 2);

  CHECK_CLOSE_ABS(vpab.getAxis().getMin(), 0.0, 1e-15);

  CHECK_CLOSE_ABS(vpab.getAxis().getMax(), 10.0, 1e-15);

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
