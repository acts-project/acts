// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/SphericalSpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"

#include <cmath>
#include <numbers>

namespace Acts::Test {

namespace {

// Minimal valid configuration for the (phi, eta, r) grid. bFieldInZ = 0 makes
// the constructor use maxPhiBins directly, so we do not need a physical helix
// sizing for the phi axis in this unit test. Four eta bins: [-3,-1], [-1,0],
// [0,1], [1,3]; two r bins: [0,100], [100,200].
Acts::Experimental::SphericalSpacePointGrid::Config makeConfig() {
  Acts::Experimental::SphericalSpacePointGrid::Config cfg;
  cfg.bFieldInZ = 0;
  cfg.maxPhiBins = 20;
  cfg.phiMin = -std::numbers::pi_v<float>;
  cfg.phiMax = std::numbers::pi_v<float>;
  cfg.rMin = 0.f;
  cfg.rMax = 200.f;
  cfg.etaMin = -3.f;
  cfg.etaMax = 3.f;
  cfg.etaBinEdges = {-3.f, -1.f, 0.f, 1.f, 3.f};
  cfg.rBinEdges = {0.f, 100.f, 200.f};
  cfg.bottomBinFinder.emplace(1, 1, 0);
  cfg.topBinFinder.emplace(1, 1, 0);
  return cfg;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SphericalSpacePointGridTests)

// The grid bins on cot(theta) = z / r = sinh(eta) and stores sinh(etaBinEdges)
// on the axis, so a space point falls into exactly the eta bin it belongs to.
// We exercise that by looking up by cot(theta) = sinh(eta).
BOOST_AUTO_TEST_CASE(BinsByCotThetaMatchEta) {
  Acts::Experimental::SphericalSpacePointGrid grid(makeConfig());

  BOOST_CHECK_GT(grid.numberOfBins(), 0u);

  const float phi = 0.f;
  const float r = 150.f;
  auto binAtEta = [&](float eta) {
    return grid.binIndex(phi, std::sinh(eta), r);
  };

  BOOST_REQUIRE(binAtEta(0.5f).has_value());

  // Two points in the same eta bin ([0,1]) map to the same global bin.
  BOOST_CHECK(binAtEta(0.5f) == binAtEta(0.7f));
  // Points in different eta bins map to different global bins.
  BOOST_CHECK(binAtEta(0.5f) != binAtEta(1.5f));   // [0,1] vs [1,3]
  BOOST_CHECK(binAtEta(0.5f) != binAtEta(-0.5f));  // [0,1] vs [-1,0]
  // The sinh mapping places the bin boundary exactly at the eta edge (eta = 1).
  BOOST_CHECK(binAtEta(0.99f) != binAtEta(1.01f));
  BOOST_CHECK(binAtEta(1.01f) == binAtEta(2.0f));  // both in [1,3]

  // The r axis separates points at the same (phi, eta) but different r.
  BOOST_CHECK(grid.binIndex(phi, std::sinh(0.5f), 50.f) !=
              grid.binIndex(phi, std::sinh(0.5f), 150.f));
}

// Points outside the eta acceptance (|eta| > etaMax) fall outside the grid and
// yield no bin; inserting such a point is a no-op for the counter.
BOOST_AUTO_TEST_CASE(OutOfAcceptanceHasNoBin) {
  Acts::Experimental::SphericalSpacePointGrid grid(makeConfig());

  const float phi = 0.f;
  const float r = 150.f;
  // eta = 5 is beyond etaMax = 3.
  BOOST_CHECK(!grid.binIndex(phi, std::sinh(5.f), r).has_value());

  const std::size_t before = grid.numberOfSpacePoints();
  const auto inserted = grid.insert(99, phi, std::sinh(5.f), r);
  BOOST_CHECK(!inserted.has_value());
  BOOST_CHECK_EQUAL(grid.numberOfSpacePoints(), before);
}

// When etaBinEdges is empty the grid builds equidistant eta bins from
// etaMin/etaMax with deltaEtaMax width (the default production path). It must
// still bin in eta.
BOOST_AUTO_TEST_CASE(UniformEtaFallback) {
  auto cfg = makeConfig();
  cfg.etaBinEdges.clear();
  cfg.deltaEtaMax = 0.5f;  // etaMin=-3 .. etaMax=3 -> 12 bins
  Acts::Experimental::SphericalSpacePointGrid grid(cfg);

  BOOST_CHECK_GT(grid.numberOfBins(), 0u);
  const float phi = 0.f;
  const float r = 150.f;
  BOOST_CHECK(grid.binIndex(phi, std::sinh(-2.f), r) !=
              grid.binIndex(phi, std::sinh(2.f), r));
}

// Inserted space points land in the expected bins: two points in the same eta
// bin share a bin, a third in a different eta bin is separate.
BOOST_AUTO_TEST_CASE(InsertPlacesInExpectedBins) {
  Acts::Experimental::SphericalSpacePointGrid grid(makeConfig());

  const float phi = 0.f;
  const float r = 150.f;
  grid.insert(0, phi, std::sinh(0.5f), r);  // eta bin [0,1]
  grid.insert(1, phi, std::sinh(0.7f), r);  // same bin
  grid.insert(2, phi, std::sinh(1.5f), r);  // eta bin [1,3]

  BOOST_CHECK_EQUAL(grid.numberOfSpacePoints(), 3u);

  const auto binLow = grid.binIndex(phi, std::sinh(0.5f), r);
  const auto binHigh = grid.binIndex(phi, std::sinh(1.5f), r);
  BOOST_REQUIRE(binLow.has_value());
  BOOST_REQUIRE(binHigh.has_value());
  BOOST_CHECK(binLow != binHigh);
  BOOST_CHECK_EQUAL(grid.at(*binLow).size(), 2u);   // indices 0 and 1
  BOOST_CHECK_EQUAL(grid.at(*binHigh).size(), 1u);  // index 2
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
