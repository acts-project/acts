// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Vertexing/DummyVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <system_error>
#include <vector>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation, i.e.
/// tests if reordering tracks returns the same result
///
BOOST_AUTO_TEST_CASE(track_density_finder_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);
  using Finder =
      TrackDensityVertexFinder<DummyVertexFitter<>,
                               GaussianTrackDensity<BoundTrackParameters>>;
  Finder finder;
  Finder::State state;

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  auto params1a = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1a, 0), mom1a,
                                               mom1a.norm(), 1, covMat)
                      .value();
  auto params1b = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1b, 0), mom1b,
                                               mom1b.norm(), -1, covMat)
                      .value();
  auto params1c = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1c, 0), mom1c,
                                               mom1c.norm(), -1, covMat)
                      .value();

  // Vectors of track parameters in different orders
  std::vector<const BoundTrackParameters*> vec1 = {&params1a, &params1b,
                                                   &params1c};
  std::vector<const BoundTrackParameters*> vec2 = {&params1c, &params1a,
                                                   &params1b};

  auto res1 = finder.find(vec1, vertexingOptions, state);
  auto res2 = finder.find(vec2, vertexingOptions, state);

  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }

  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }

  if (res1.ok() and res2.ok()) {
    BOOST_CHECK(!(*res1).empty());
    BOOST_CHECK(!(*res2).empty());
    Vector3 result1 = (*res1).back().position();
    Vector3 result2 = (*res2).back().position();
    BOOST_CHECK(result1 == result2);
  }
}

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation
///
BOOST_AUTO_TEST_CASE(track_density_finder_constr_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Finder options
  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  // Create constraint for seed finding
  Vector3 constraintPos{1.7_mm, 1.3_mm, -6_mm};
  SymMatrix3 constrCov = ActsSymMatrix<3>::Identity();

  Vertex<BoundTrackParameters> vertexConstraint(constraintPos);
  vertexConstraint.setCovariance(constrCov);

  vertexingOptions.vertexConstraint = vertexConstraint;
  using Finder =
      TrackDensityVertexFinder<DummyVertexFitter<>,
                               GaussianTrackDensity<BoundTrackParameters>>;
  Finder finder;
  Finder::State state;

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  auto params1a = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1a, 0), mom1a,
                                               mom1a.norm(), 1, covMat)
                      .value();
  auto params1b = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1b, 0), mom1b,
                                               mom1b.norm(), -1, covMat)
                      .value();
  auto params1c = BoundTrackParameters::create(perigeeSurface, geoContext,
                                               makeVector4(pos1c, 0), mom1c,
                                               mom1c.norm(), -1, covMat)
                      .value();

  // Vector of track parameters
  std::vector<const BoundTrackParameters*> vec1 = {&params1a, &params1b,
                                                   &params1c};

  auto res = finder.find(vec1, vertexingOptions, state);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3 result = (*res).back().position();

    BOOST_CHECK(result[eX] == constraintPos[eX]);
    BOOST_CHECK(result[eY] == constraintPos[eY]);
    CHECK_CLOSE_ABS(result[eZ], expectedZResult, 0.001_mm);
  }
}

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
BOOST_AUTO_TEST_CASE(track_density_finder_random_test) {
  Covariance covMat = Covariance::Identity();

  // Perigee surface for track parameters
  Vector3 pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);
  using Finder =
      TrackDensityVertexFinder<DummyVertexFitter<>,
                               GaussianTrackDensity<BoundTrackParameters>>;
  Finder finder;
  Finder::State state;

  int mySeed = 31415;
  std::mt19937 gen(mySeed);
  unsigned int nTracks = 200;

  std::vector<BoundTrackParameters> trackVec;
  trackVec.reserve(nTracks);

  // Create nTracks tracks for test case
  for (unsigned int i = 0; i < nTracks; i++) {
    // The position of the particle
    Vector3 pos(xdist(gen), ydist(gen), 0);

    // Create momentum and charge of track
    double pt = pTDist(gen);
    double phi = phiDist(gen);
    double eta = etaDist(gen);
    double charge = etaDist(gen) > 0 ? 1 : -1;

    // project the position on the surface
    Vector3 direction = makeDirectionUnitFromPhiEta(phi, eta);
    auto intersection = perigeeSurface->intersect(geoContext, pos, direction);
    pos = intersection.intersection.position;

    // Produce most of the tracks at near z1 position,
    // some near z2. Highest track density then expected at z1
    pos[eZ] = ((i % 4) == 0) ? z2dist(gen) : z1dist(gen);

    trackVec.push_back(BoundTrackParameters::create(
                           perigeeSurface, geoContext, makeVector4(pos, 0),
                           direction, pt, charge, covMat)
                           .value());
  }

  std::vector<const BoundTrackParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res3 = finder.find(trackPtrVec, vertexingOptions, state);
  if (!res3.ok()) {
    std::cout << res3.error().message() << std::endl;
  }

  if (res3.ok()) {
    BOOST_CHECK(!(*res3).empty());
    Vector3 result = (*res3).back().position();
    CHECK_CLOSE_ABS(result[eZ], zVertexPos, 1_mm);
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundTrackParameters& params) : m_parameters(params) {}

  const BoundTrackParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundTrackParameters m_parameters;
};

///
/// @brief Unit test for TrackDensityVertexFinder with user-defined input track
/// type with same values as in other tests
///
BOOST_AUTO_TEST_CASE(track_density_finder_usertrack_test) {
  // Define some track parameter properties
  Vector3 pos0{0, 0, 0};
  Vector3 pos1a{1.86052_mm, -1.24035_mm, -10_mm};
  Vector3 mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3 pos1b{-1.24035_mm, 1.86052_mm, -3_mm};
  Vector3 mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3 pos1c{1.69457_mm, -0.50837_mm, -7_mm};
  Vector3 mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Finder options
  VertexingOptions<InputTrack> vertexingOptions(geoContext, magFieldContext);

  // Create constraint for seed finding
  Vector3 constraintPos{1.7_mm, 1.3_mm, -6_mm};
  SymMatrix3 constrCov = SymMatrix3::Identity();

  Vertex<InputTrack> vertexConstraint(constraintPos);
  vertexConstraint.setCovariance(constrCov);

  vertexingOptions.vertexConstraint = vertexConstraint;

  std::function<BoundTrackParameters(InputTrack)> extractParameters =
      [](const InputTrack& params) { return params.parameters(); };

  using Finder = TrackDensityVertexFinder<DummyVertexFitter<InputTrack>,
                                          GaussianTrackDensity<InputTrack>>;

  Finder finder(extractParameters);
  Finder::State state;

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  InputTrack params1a(BoundTrackParameters::create(perigeeSurface, geoContext,
                                                   makeVector4(pos1a, 0), mom1a,
                                                   mom1a.norm(), 1, covMat)
                          .value());
  InputTrack params1b(BoundTrackParameters::create(perigeeSurface, geoContext,
                                                   makeVector4(pos1b, 0), mom1b,
                                                   mom1b.norm(), -1, covMat)
                          .value());
  InputTrack params1c(BoundTrackParameters::create(perigeeSurface, geoContext,
                                                   makeVector4(pos1c, 0), mom1c,
                                                   mom1c.norm(), -1, covMat)
                          .value());

  // Vector of track parameters
  std::vector<const InputTrack*> vec1 = {&params1a, &params1b, &params1c};

  auto res = finder.find(vec1, vertexingOptions, state);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3 result = (*res).back().position();

    BOOST_CHECK(result[eX] == constraintPos[eX]);
    BOOST_CHECK(result[eY] == constraintPos[eY]);
    CHECK_CLOSE_ABS(result[eZ], expectedZResult, 0.001_mm);
  }
}

}  // namespace Test
}  // namespace Acts
