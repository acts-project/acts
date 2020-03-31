// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

const double zVertexPos = 12.;
// x position
std::normal_distribution<double> xdist(1_mm, 0.1_mm);
// y position
std::normal_distribution<double> ydist(-0.7_mm, 0.1_mm);
// z1 position
std::normal_distribution<double> z1dist(zVertexPos * 1_mm, 1_mm);
// z2 position
std::normal_distribution<double> z2dist(-3_mm, 0.5_mm);
// Track pT distribution
std::uniform_real_distribution<double> pTDist(0.1_GeV, 100_GeV);
// Track phi distribution
std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
// Track eta distribution
std::uniform_real_distribution<double> etaDist(-4., 4.);

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation
///
BOOST_AUTO_TEST_CASE(grid_density_vertex_finder_test) {
  const int mainGridSize = 3000;
  const int trkGridSize = 35;

  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3D pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundParameters> vertexingOptions(geoContext,
                                                     magFieldContext);

  GridDensityVertexFinder<mainGridSize, trkGridSize>::Config cfg;
  GridDensityVertexFinder<mainGridSize, trkGridSize> finder(cfg);

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3D pos(xdist(gen), ydist(gen), 0);
    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    if ((i % 4) == 0) {
      pos[eZ] = z2dist(gen);
    } else {
      pos[eZ] = z1dist(gen);
    }

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    Vector3D mom(pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta));
    double charge = etaDist(gen) > 0 ? 1 : -1;

    trackVec.push_back(BoundParameters(geoContext, covMat, pos, mom, charge, 0,
                                       perigeeSurface));
  }

  std::vector<const BoundParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res3 = finder.find(trackPtrVec, vertexingOptions);
  if (!res3.ok()) {
    std::cout << res3.error().message() << std::endl;
  }

  if (res3.ok()) {
    BOOST_CHECK(!(*res3).empty());
    Vector3D result = (*res3).back().position();
    std::cout << "result: " << result << std::endl;
    CHECK_CLOSE_ABS(result[eZ], zVertexPos, 1_mm);
  }
}

}  // namespace Test
}  // namespace Acts
