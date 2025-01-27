// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinUtility.hpp"

using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(ProtoBinningPlaceHolderEquidistant) {
  // A valid placeholder binning
  auto peq =
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, 5u);
  BOOST_CHECK_EQUAL(peq.bins(), 5u);
}

BOOST_AUTO_TEST_CASE(ProtoBinningEquidistant) {
  // An invalid binning, 0 bins given
  BOOST_CHECK_THROW(
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, 15., 20.,
                   0),
      std::invalid_argument);

  // Another invalid binning, min/max swapped
  BOOST_CHECK_THROW(
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, 150., 20.,
                   10),
      std::invalid_argument);

  // A valid binning
  auto eq = ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, 0.,
                         10., 5u);

  std::vector<Acts::ActsScalar> reference = {0., 2., 4., 6., 8., 10.};
  BOOST_CHECK_EQUAL(eq.bins(), 5u);
  BOOST_CHECK_EQUAL(eq.binValue, Acts::binX);
  BOOST_CHECK(eq.axisType == Acts::detail::AxisType::Equidistant);
  BOOST_CHECK(eq.boundaryType == Acts::detail::AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL_COLLECTIONS(eq.edges.begin(), eq.edges.end(),
                                reference.begin(), reference.end());
}

BOOST_AUTO_TEST_CASE(ProtoBinningVariable) {
  // An invalid binning, edge size < 2u
  std::vector<Acts::ActsScalar> iedges = {12.};
  BOOST_CHECK_THROW(
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, iedges),
      std::invalid_argument);

  // A valid binning
  std::vector<Acts::ActsScalar> varEdges = {0., 12., 13., 15., 20.};
  auto var =
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, varEdges);

  BOOST_CHECK_EQUAL(var.bins(), 4u);
  BOOST_CHECK_EQUAL(var.binValue, Acts::binX);
  BOOST_CHECK(var.axisType == Acts::detail::AxisType::Variable);
  BOOST_CHECK(var.boundaryType == Acts::detail::AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL_COLLECTIONS(var.edges.begin(), var.edges.end(),
                                varEdges.begin(), varEdges.end());
}

BOOST_AUTO_TEST_CASE(BinningDescriptionFromAndToBinUtility) {
  // A valid binning
  Acts::BinUtility bUtility(5u, 0., 10., Acts::open, Acts::binR);
  std::vector<float> edges = {-M_PI, 0.1, M_PI};
  bUtility += Acts::BinUtility(edges, Acts::closed, Acts::binPhi);

  auto bDescription = BinningDescription::fromBinUtility(bUtility);

  BOOST_CHECK_EQUAL(bDescription.binning.size(), 2u);

  // Test the first entry
  BOOST_CHECK_EQUAL(bDescription.binning[0].bins(), 5u);
  BOOST_CHECK_EQUAL(bDescription.binning[0].binValue, Acts::binR);
  BOOST_CHECK(bDescription.binning[0].axisType ==
              Acts::detail::AxisType::Equidistant);
  BOOST_CHECK(bDescription.binning[0].boundaryType ==
              Acts::detail::AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL(bDescription.binning[0].edges.size(), 6u);

  // Check the second entry
  BOOST_CHECK_EQUAL(bDescription.binning[1].bins(), 2u);
  BOOST_CHECK_EQUAL(bDescription.binning[1].binValue, Acts::binPhi);
  BOOST_CHECK(bDescription.binning[1].axisType ==
              Acts::detail::AxisType::Variable);
  BOOST_CHECK(bDescription.binning[1].boundaryType ==
              Acts::detail::AxisBoundaryType::Closed);
  BOOST_CHECK_EQUAL(bDescription.binning[1].edges.size(), 3u);

  // Round-trip
  auto binUtility = bDescription.toBinUtility();
  BOOST_CHECK_EQUAL(binUtility.binningData().size(), 2u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[0].bins(), 5u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[1].bins(), 2u);
  BOOST_CHECK_EQUAL(binUtility.binningData()[1].boundaries().size(), 3u);
  CHECK_CLOSE_ABS(binUtility.binningData()[1].boundaries()[0], -M_PI, 1e-5);
  CHECK_CLOSE_ABS(binUtility.binningData()[1].boundaries()[1], 0.1, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
