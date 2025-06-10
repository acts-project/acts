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
  using enum Acts::AxisType;

  // Bound, equidistant axis
  Acts::ProtoAxis epab(Bound, 0.0, 1.0, 10);

  // Direct access
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

  std::string rString = "ProtoAxis: 10 bins, equidistant within [0, 1]";
  std::string oString = epab.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Test copy constructor
  Acts::ProtoAxis epabCopy(epab);
  BOOST_CHECK_EQUAL(epabCopy.isAutorange(), epab.isAutorange());
  BOOST_CHECK(epabCopy.getAxis() == epab.getAxis());

  // Test Assignment operator
  Acts::ProtoAxis epabAssign(Bound, 0.0, 50.0, 20);
  epabAssign = epab;
  BOOST_CHECK_EQUAL(epabAssign.isAutorange(), epab.isAutorange());
  BOOST_CHECK(epabAssign.getAxis() == epab.getAxis());

  // Create a grid from a single proto axis
  auto grid1D = Acts::makeGrid<double>(epab);
  BOOST_CHECK(grid1D != nullptr);
  BOOST_CHECK_EQUAL(grid1D->axes().size(), 1);
  auto axis1D =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Bound>*>(
          grid1D->axes().front());
  BOOST_CHECK(axis1D != nullptr);

  // Open, equidistant axis
  Acts::ProtoAxis epao(Open, 0., 2.0, 10.);
  BOOST_CHECK_EQUAL(epao.getAxis().getBoundaryType(), Open);

  // Create a 2D grid from a two proto axes
  auto grid2D = Acts::makeGrid<double>(epab, epao);
  BOOST_CHECK(grid2D != nullptr);
  auto grid2Daxes = grid2D->axes();
  BOOST_CHECK_EQUAL(grid2Daxes.size(), 2);
  auto axis2D1 =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Bound>*>(
          grid2Daxes[0]);
  BOOST_CHECK(axis2D1 != nullptr);
  auto axis2D2 =
      dynamic_cast<const Acts::Axis<Acts::AxisType::Equidistant, Open>*>(
          grid2Daxes[1]);
  BOOST_CHECK(axis2D2 != nullptr);

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(Bound, 1.0, 0.0, 10),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(Bound, 0.0, 1.0, 0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AutorangeProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisType;

  // Bound, equidistant axis, autorange
  Acts::ProtoAxis epa(Bound, 10);

  // Direct access
  BOOST_CHECK(epa.isAutorange());

  // Access via IAxis
  BOOST_CHECK(epa.getAxis().isEquidistant());

  BOOST_CHECK(!epa.getAxis().isVariable());

  BOOST_CHECK_EQUAL(epa.getAxis().getType(), Equidistant);

  BOOST_CHECK_EQUAL(epa.getAxis().getBoundaryType(), Bound);

  BOOST_CHECK_EQUAL(epa.getAxis().getNBins(), 10);

  std::string rString =
      "ProtoAxis: 10 bins, equidistant within automatic range";
  std::string oString = epa.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Test copy constructor
  Acts::ProtoAxis epaCopy(epa);
  BOOST_CHECK_EQUAL(epaCopy.isAutorange(), epa.isAutorange());
  BOOST_CHECK(epaCopy.getAxis() == epa.getAxis());

  // Test Assignment operator
  Acts::ProtoAxis epaAssign(Bound, 0.0, 50.0, 20);
  epaAssign = epa;
  BOOST_CHECK_EQUAL(epaAssign.isAutorange(), epa.isAutorange());
  BOOST_CHECK(epaAssign.getAxis() == epa.getAxis());

  // Invalid 1D grid construction with autorange axis
  BOOST_CHECK_THROW(Acts::makeGrid<double>(epa), std::invalid_argument);

  // Invalid 2D grid construction with autorange axis
  Acts::ProtoAxis epao(Open, 0., 2.0, 10.);
  BOOST_CHECK_THROW(Acts::makeGrid<double>(epao, epa), std::invalid_argument);
  BOOST_CHECK_THROW(Acts::makeGrid<double>(epa, epao), std::invalid_argument);

  // Set the range now
  epa.setRange(0.0, 20.0);
  BOOST_CHECK(!epa.isAutorange());

  // 1D Grid consstruction works now
  BOOST_CHECK_NO_THROW(Acts::makeGrid<double>(epa));

  // 2D Grid consstruction works now
  BOOST_CHECK_NO_THROW(Acts::makeGrid<double>(epa, epao));

  // Invalid setRange, min > max
  BOOST_CHECK_THROW(epa.setRange(20.0, 0.0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VariableProtoAxis) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisType;

  // Bound, variable axis
  Acts::ProtoAxis vpab(Bound, {0.0, 1.0, 10});

  // Direct access
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

  std::string rString = "ProtoAxis: 2 bins, variable within [0, 10]";
  std::string oString = vpab.toString();
  BOOST_CHECK_EQUAL(rString, oString);

  // Test copy constructor
  Acts::ProtoAxis vpabCopy(vpab);
  BOOST_CHECK_EQUAL(vpabCopy.isAutorange(), vpab.isAutorange());
  BOOST_CHECK(vpabCopy.getAxis() == vpab.getAxis());

  // Test Assignment operator
  Acts::ProtoAxis vpabAssign(Bound, 0.0, 50.0, 20);
  vpabAssign = vpab;
  BOOST_CHECK_EQUAL(vpabAssign.isAutorange(), vpab.isAutorange());
  BOOST_CHECK(vpabAssign.getAxis() == vpab.getAxis());

  // Set Range for variable
  vpab.setRange(0.5, 9.5);
  // Bins stay the same. min/max update
  BOOST_CHECK_EQUAL(vpab.getAxis().getNBins(), 2);
  CHECK_CLOSE_ABS(vpab.getAxis().getMin(), 0.5, 1e-15);
  CHECK_CLOSE_ABS(vpab.getAxis().getMax(), 9.5, 1e-15);

  // Invalid constructor, min > max
  BOOST_CHECK_THROW(Acts::ProtoAxis(Bound, std::vector<double>{2.}),
                    std::invalid_argument);

  // Invalid constructor, nbins < 1
  BOOST_CHECK_THROW(Acts::ProtoAxis(Bound, {3., 2., 1}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DirectedProtoAxisOfAllSorts) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;
  using enum Acts::AxisType;

  // Bound, variable axis - directed at x
  Acts::DirectedProtoAxis vpab(AxisX, Bound, {0.0, 1.0, 10});
  BOOST_CHECK(vpab.getAxisDirection() == AxisX);
  std::string rvString =
      "DirectedProtoAxis: 2 bins in AxisX, variable within [0, 10]";
  std::string ovString = vpab.toString();
  BOOST_CHECK_EQUAL(rvString, ovString);

  // Bound, equidistant axis, autorange - directed at y
  Acts::DirectedProtoAxis epa(AxisY, Bound, 10);
  BOOST_CHECK(epa.getAxisDirection() == AxisY);
  BOOST_CHECK(epa.getAxis().isEquidistant());
  BOOST_CHECK(epa.isAutorange());
  std::string reString =
      "DirectedProtoAxis: 10 bins in AxisY, equidistant within automatic range";
  std::string oeString = epa.toString();
  BOOST_CHECK_EQUAL(reString, oeString);

  // Bound, equidistant axis - directed at z
  Acts::DirectedProtoAxis epab(AxisZ, Bound, 0.0, 1.0, 10);
  BOOST_CHECK(epab.getAxisDirection() == AxisZ);
  BOOST_CHECK(epab.getAxis().isEquidistant());
  BOOST_CHECK(!epab.isAutorange());
  std::string rString =
      "DirectedProtoAxis: 10 bins in AxisZ, equidistant within [0, 1]";
  std::string oString = epab.toString();
  BOOST_CHECK_EQUAL(rString, oString);
}

BOOST_AUTO_TEST_SUITE_END()
