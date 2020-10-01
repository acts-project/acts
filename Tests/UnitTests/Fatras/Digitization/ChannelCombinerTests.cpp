// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Digitization/ChannelCombiner.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

#include <Acts/Tests/CommonHelpers/FloatComparisons.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(ChannelCombiner1D) {
  // The channel combiner in 1D
  ChannelCombiner cc;

  // Some cells
  Cell cell0(5, 10.5);
  Cell cell1(6, 11.5);
  Cell cell2(7, 12.5);

  // Digital clustering test
  std::vector<Channel<double, Acts::eBoundLoc0>> channelsD = {
      {{cell0}, 1.}, {{cell1}, 1.}, {{cell2}, 1.}};

  auto digitalCluster = cc.combine<double, Acts::eBoundLoc0>(channelsD);

  CHECK_CLOSE_ABS(
      digitalCluster.parameterSet.template getParameter<Acts::eBoundLoc0>(),
      11.5, Acts::s_epsilon);

  auto loc0Index =
      digitalCluster.parameterSet.template getIndex<Acts::eBoundLoc0>();
  CHECK_CLOSE_ABS(digitalCluster.clusterSize[loc0Index], 3, Acts::s_epsilon);

  // Analog clustering test
  std::vector<Channel<double, Acts::eBoundLoc0>> channelsA = {
      {{cell0}, 0.5}, {{cell1}, 1.}, {{cell2}, 0.25}};

  auto analogCluster = cc.combine<double, Acts::eBoundLoc0>(channelsA);

  CHECK_CLOSE_ABS(
      analogCluster.parameterSet.template getParameter<Acts::eBoundLoc0>(),
      11.35714, 0.001);

  CHECK_CLOSE_ABS(analogCluster.clusterSize[loc0Index], 3, Acts::s_epsilon);
}

BOOST_AUTO_TEST_CASE(ChannelCombiner2D) {
  // The channel combiner in 1D
  ChannelCombiner cc;

  // Some cells
  Cell cell00(5, 10.5);
  Cell cell01(6, 11.5);
  Cell cell02(7, 12.5);

  Cell cell10(2, -3.5);
  Cell cell11(2, -3.6);
  Cell cell12(3, -3.7);

  // Digital clustering test (analog is test in 1D already)
  std::vector<Channel<double, Acts::eBoundLoc0, Acts::eBoundLoc1>> channelsD = {
      {{cell00, cell10}, 1.}, {{cell01, cell11}, 1.}, {{cell02, cell12}, 1.}};

  auto digitalCluster =
      cc.combine<double, Acts::eBoundLoc0, Acts::eBoundLoc1>(channelsD);
  CHECK_CLOSE_ABS(
      digitalCluster.parameterSet.template getParameter<Acts::eBoundLoc0>(),
      11.5, Acts::s_epsilon);
  CHECK_CLOSE_ABS(
      digitalCluster.parameterSet.template getParameter<Acts::eBoundLoc1>(),
      -3.6, Acts::s_epsilon);

  auto loc0Index =
      digitalCluster.parameterSet.template getIndex<Acts::eBoundLoc0>();
  auto loc1Index =
      digitalCluster.parameterSet.template getIndex<Acts::eBoundLoc1>();

  CHECK_CLOSE_ABS(digitalCluster.clusterSize[loc0Index], 3, Acts::s_epsilon);
  CHECK_CLOSE_ABS(digitalCluster.clusterSize[loc1Index], 2, Acts::s_epsilon);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
