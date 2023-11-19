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

using namespace Acts::Experimental;

BOOST_AUTO_TEST_SUITE(Detector)

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
  BOOST_CHECK_THROW(
      ProtoBinning(Acts::binX, Acts::detail::AxisBoundaryType::Bound, {12.}),
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

BOOST_AUTO_TEST_SUITE_END()
