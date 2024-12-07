// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <numbers>

using namespace Acts;
using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(ProtoBinningPlaceHolderEquidistant) {
  // A valid placeholder binning
  auto peq = ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound, 5u);
  BOOST_CHECK_EQUAL(peq.bins(), 5u);
}

BOOST_AUTO_TEST_CASE(ProtoBinningEquidistant) {
  // An invalid binning, 0 bins given
  BOOST_CHECK_THROW(
      ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound, 15., 20., 0),
      std::invalid_argument);

  // Another invalid binning, min/max swapped
  BOOST_CHECK_THROW(ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound,
                                 150., 20., 10),
                    std::invalid_argument);

  // A valid binning
  auto eq =
      ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound, 0., 10., 5u);

  std::vector<double> reference = {0., 2., 4., 6., 8., 10.};
  BOOST_CHECK_EQUAL(eq.bins(), 5u);
  BOOST_CHECK_EQUAL(eq.axisDirection, AxisDirection::AxisX);
  BOOST_CHECK(eq.axisType == AxisType::Equidistant);
  BOOST_CHECK(eq.axisBoundaryType == AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL_COLLECTIONS(eq.edges.begin(), eq.edges.end(),
                                reference.begin(), reference.end());
}

BOOST_AUTO_TEST_CASE(ProtoBinningVariable) {
  // An invalid binning, edge size < 2u
  std::vector<double> iedges = {12.};
  BOOST_CHECK_THROW(
      ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound, iedges),
      std::invalid_argument);

  // A valid binning
  std::vector<double> varEdges = {0., 12., 13., 15., 20.};
  auto var =
      ProtoBinning(AxisDirection::AxisX, AxisBoundaryType::Bound, varEdges);

  BOOST_CHECK_EQUAL(var.bins(), 4u);
  BOOST_CHECK_EQUAL(var.axisDirection, AxisDirection::AxisX);
  BOOST_CHECK(var.axisType == AxisType::Variable);
  BOOST_CHECK(var.axisBoundaryType == AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL_COLLECTIONS(var.edges.begin(), var.edges.end(),
                                varEdges.begin(), varEdges.end());
}

BOOST_AUTO_TEST_CASE(BinningDescriptionFromAndToBinUtility) {
  // A valid binning
  Acts::BinUtility bUtility(5u, 0., 10., Acts::AxisBoundaryType::Bound,
                            AxisDirection::AxisR);
  std::vector<float> edges = {-std::numbers::pi, 0.1, std::numbers::pi};
  bUtility += Acts::BinUtility(edges, Acts::AxisBoundaryType::Closed,
                               AxisDirection::AxisPhi);

  auto bDescription = BinningDescription::fromBinUtility(bUtility);

  BOOST_CHECK_EQUAL(bDescription.binning.size(), 2u);

  // Test the first entry
  BOOST_CHECK_EQUAL(bDescription.binning[0].bins(), 5u);
  BOOST_CHECK_EQUAL(bDescription.binning[0].axisDirection,
                    AxisDirection::AxisR);
  BOOST_CHECK(bDescription.binning[0].axisType == AxisType::Equidistant);
  BOOST_CHECK(bDescription.binning[0].axisBoundaryType ==
              AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL(bDescription.binning[0].edges.size(), 6u);

  // Check the second entry
  BOOST_CHECK_EQUAL(bDescription.binning[1].bins(), 2u);
  BOOST_CHECK_EQUAL(bDescription.binning[1].axisDirection,
                    AxisDirection::AxisPhi);
  BOOST_CHECK(bDescription.binning[1].axisType == AxisType::Variable);
  BOOST_CHECK(bDescription.binning[1].axisBoundaryType ==
              AxisBoundaryType::Closed);
  BOOST_CHECK_EQUAL(bDescription.binning[1].edges.size(), 3u);

  // Round-trip
  auto binUtility = bDescription.toBinUtility();
  BOOST_CHECK_EQUAL(binUtility.binningData().size(), 2u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[0].bins(), 5u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[1].bins(), 2u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[1].boundaries().size(), 3u);
  CHECK_CLOSE_ABS(binUtility.binningData()[1].boundaries()[0],
                  -std::numbers::pi, 1e-5);
  CHECK_CLOSE_ABS(binUtility.binningData()[1].boundaries()[1], 0.1, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
