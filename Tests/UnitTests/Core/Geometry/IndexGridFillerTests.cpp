// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/IndexGrid.hpp"

#include <vector>

using namespace Acts;

auto tContext = GeometryContext::dangerouslyDefaultConstruct();
Logging::Level logLevel = Logging::VERBOSE;

namespace {

/// Helper method to count how many bins are not empty
template <typename indexed_surface_grid>
std::size_t countBins(const indexed_surface_grid& isGrid) {
  std::size_t nonEmptyBins = 0u;
  for (std::size_t igb = 0u; igb < isGrid.grid.size(); ++igb) {
    const auto& gb = isGrid.grid.at(igb);
    if (!gb.empty()) {
      ++nonEmptyBins;
    }
  }
  return nonEmptyBins;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(DetectorSuite)

BOOST_AUTO_TEST_CASE(BinSequence) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Pre-Test", logLevel));
  ACTS_INFO("Testing bin sequence generators.");

  // Test standard bound local bin sequence
  auto seq48e0b10B = binSequence({4u, 8u}, 0u, 10u, AxisBoundaryType::Bound);
  std::vector<std::size_t> reference = {4u, 5u, 6u, 7u, 8u};
  BOOST_CHECK(seq48e0b10B == reference);

  // Test bound local bin sequence with expansion 1u
  auto seq48e1b10B = binSequence({4u, 8u}, 1u, 10u, AxisBoundaryType::Bound);
  reference = {3u, 4u, 5u, 6u, 7u, 8u, 9u};
  BOOST_CHECK(seq48e1b10B == reference);

  // Test bound local bin sequence with expansion 3u - clipped to max bin 10u
  auto seq48e3b10B = binSequence({4u, 8u}, 3u, 10u, AxisBoundaryType::Bound);
  reference = {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u};
  BOOST_CHECK(seq48e3b10B == reference);

  // Test open bin sequence with overflow filling
  auto seq48e3b10O = binSequence({4u, 8u}, 3u, 10u, AxisBoundaryType::Open);
  reference = {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u};
  BOOST_CHECK(seq48e3b10O == reference);

  // Test standard closed local bins
  auto seq48e0b10C = binSequence({4u, 8u}, 0u, 20u, AxisBoundaryType::Closed);
  reference = {4u, 5u, 6u, 7u, 8u};
  BOOST_CHECK(seq48e0b10C == reference);

  // Test closed local bins with expansion
  auto seq48e1b10C = binSequence({4u, 8u}, 1u, 20u, AxisBoundaryType::Closed);
  reference = {3u, 4u, 5u, 6u, 7u, 8u, 9u};
  BOOST_CHECK(seq48e1b10C == reference);

  // Test closed local bins with expansion over bin boundary
  auto seq1029e1b20C =
      binSequence({19u, 20u}, 1u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 18u, 19u, 20u};
  BOOST_CHECK(seq1029e1b20C == reference);

  // Test closed local bins with bin boundary jump
  auto seq218e0b20C = binSequence({2u, 18u}, 0u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 2u, 18u, 19u, 20u};
  BOOST_CHECK(seq218e0b20C == reference);

  // Test closed local bins with bin boundary jump with extension
  auto seq218e2b20C = binSequence({2u, 18u}, 2u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 2u, 3u, 4u, 16u, 17u, 18u, 19u, 20u};
  BOOST_CHECK(seq218e2b20C == reference);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
