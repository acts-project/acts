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

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/DummyVertexFitter.hpp"
#include "Acts/Vertexing/GaussianTrackDensity.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation, i.e.
/// tests if reordering tracks returns the same result
///
BOOST_AUTO_TEST_CASE(track_density_finder_test) {
  // Define some track parameter properties
  Vector3D pos0{0, 0, 0};
  Vector3D pos1a{2_mm, 1_mm, -10_mm};
  Vector3D mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3D pos1b{1_mm, 2_mm, -3_mm};
  Vector3D mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3D pos1c{1.2_mm, 1.3_mm, -7_mm};
  Vector3D mom1c{300_MeV, 1000_MeV, 100_MeV};

  VertexFinderOptions<BoundParameters> vFinderOptions(tgContext, mfContext);

  TrackDensityVertexFinder<DummyVertexFitter<>, GaussianTrackDensity> finder;

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  BoundParameters params1a(tgContext, covMat, pos1a, mom1a, 1, 0,
                           perigeeSurface);
  BoundParameters params1b(tgContext, covMat, pos1b, mom1b, -1, 0,
                           perigeeSurface);
  BoundParameters params1c(tgContext, covMat, pos1c, mom1c, -1, 0,
                           perigeeSurface);

  // Vectors of track parameters in different orders
  std::vector<const BoundParameters*> vec1 = {&params1a, &params1b, &params1c};
  std::vector<const BoundParameters*> vec2 = {&params1c, &params1a, &params1b};

  auto res1 = finder.find(vec1, vFinderOptions);
  auto res2 = finder.find(vec2, vFinderOptions);

  if (!res1.ok()) {
    std::cout << res1.error().message() << std::endl;
  }

  if (!res2.ok()) {
    std::cout << res2.error().message() << std::endl;
  }

  if (res1.ok() and res2.ok()) {
    BOOST_CHECK(!(*res1).empty());
    BOOST_CHECK(!(*res2).empty());
    Vector3D result1 = (*res1).back().position();
    Vector3D result2 = (*res2).back().position();
    BOOST_CHECK(result1 == result2);
  }
}

///
/// @brief Unit test for TrackDensityVertexFinder using same configuration
/// and values as VertexSeedFinderTestAlg in Athena implementation
///
BOOST_AUTO_TEST_CASE(track_density_finder_constr_test) {
  // Define some track parameter properties
  Vector3D pos0{0, 0, 0};
  Vector3D pos1a{2_mm, 1_mm, -10_mm};
  Vector3D mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3D pos1b{1_mm, 2_mm, -3_mm};
  Vector3D mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3D pos1c{1.2_mm, 1.3_mm, -7_mm};
  Vector3D mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Finder options
  VertexFinderOptions<BoundParameters> vFinderOptions(tgContext, mfContext);

  // Create constraint for seed finding
  Vector3D constraintPos{1.7_mm, 1.3_mm, -6_mm};
  ActsSymMatrixD<3> constrCov = ActsSymMatrixD<3>::Identity();

  Vertex<BoundParameters> vertexConstraint(constraintPos);
  vertexConstraint.setCovariance(constrCov);

  vFinderOptions.vertexConstraint = vertexConstraint;

  TrackDensityVertexFinder<DummyVertexFitter<>, GaussianTrackDensity> finder;

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  BoundParameters params1a(tgContext, covMat, pos1a, mom1a, 1, 0,
                           perigeeSurface);
  BoundParameters params1b(tgContext, covMat, pos1b, mom1b, -1, 0,
                           perigeeSurface);
  BoundParameters params1c(tgContext, covMat, pos1c, mom1c, -1, 0,
                           perigeeSurface);

  // Vector of track parameters
  std::vector<const BoundParameters*> vec1 = {&params1a, &params1b, &params1c};

  auto res = finder.find(vec1, vFinderOptions);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3D result = (*res).back().position();

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
  Vector3D pos0{0, 0, 0};
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  VertexFinderOptions<BoundParameters> vFinderOptions(tgContext, mfContext);

  TrackDensityVertexFinder<DummyVertexFitter<>, GaussianTrackDensity> finder;

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

    trackVec.push_back(BoundParameters(tgContext, covMat, pos, mom, charge, 0,
                                       perigeeSurface));
  }

  std::vector<const BoundParameters*> trackPtrVec;
  for (const auto& trk : trackVec) {
    trackPtrVec.push_back(&trk);
  }

  auto res3 = finder.find(trackPtrVec, vFinderOptions);
  if (!res3.ok()) {
    std::cout << res3.error().message() << std::endl;
  }

  if (res3.ok()) {
    BOOST_CHECK(!(*res3).empty());
    Vector3D result = (*res3).back().position();
    CHECK_CLOSE_ABS(result[eZ], zVertexPos, 1_mm);
  }
}

// Dummy user-defined InputTrack type
struct InputTrack {
  InputTrack(const BoundParameters& params) : m_parameters(params) {}

  const BoundParameters& parameters() const { return m_parameters; }

  // store e.g. link to original objects here

 private:
  BoundParameters m_parameters;
};

///
/// @brief Unit test for TrackDensityVertexFinder with user-defined input track
/// type with same values as in other tests
///
BOOST_AUTO_TEST_CASE(track_density_finder_usertrack_test) {
  // Define some track parameter properties
  Vector3D pos0{0, 0, 0};
  Vector3D pos1a{2_mm, 1_mm, -10_mm};
  Vector3D mom1a{400_MeV, 600_MeV, 200_MeV};
  Vector3D pos1b{1_mm, 2_mm, -3_mm};
  Vector3D mom1b{600_MeV, 400_MeV, -200_MeV};
  Vector3D pos1c{1.2_mm, 1.3_mm, -7_mm};
  Vector3D mom1c{300_MeV, 1000_MeV, 100_MeV};

  // From Athena VertexSeedFinderTestAlg
  double const expectedZResult = -13.013;

  // Finder options
  VertexFinderOptions<InputTrack> vFinderOptions(tgContext, mfContext);

  // Create constraint for seed finding
  Vector3D constraintPos{1.7_mm, 1.3_mm, -6_mm};
  ActsSymMatrixD<3> constrCov = ActsSymMatrixD<3>::Identity();

  Vertex<InputTrack> vertexConstraint(constraintPos);
  vertexConstraint.setCovariance(constrCov);

  vFinderOptions.vertexConstraint = vertexConstraint;

  std::function<BoundParameters(InputTrack)> extractParameters =
      [](InputTrack params) { return params.parameters(); };

  TrackDensityVertexFinder<DummyVertexFitter<InputTrack>, GaussianTrackDensity>
      finder(extractParameters);

  // Start creating some track parameters
  Covariance covMat = Covariance::Identity();
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos0);

  // Test finder for some fixed track parameter values
  InputTrack params1a(
      BoundParameters(tgContext, covMat, pos1a, mom1a, 1, 0, perigeeSurface));
  InputTrack params1b(
      BoundParameters(tgContext, covMat, pos1b, mom1b, -1, 0, perigeeSurface));
  InputTrack params1c(
      BoundParameters(tgContext, covMat, pos1c, mom1c, -1, 0, perigeeSurface));

  // Vector of track parameters
  std::vector<const InputTrack*> vec1 = {&params1a, &params1b, &params1c};

  auto res = finder.find(vec1, vFinderOptions);

  if (!res.ok()) {
    std::cout << res.error().message() << std::endl;
  }

  if (res.ok()) {
    BOOST_CHECK(!(*res).empty());
    Vector3D result = (*res).back().position();

    BOOST_CHECK(result[eX] == constraintPos[eX]);
    BOOST_CHECK(result[eY] == constraintPos[eY]);
    CHECK_CLOSE_ABS(result[eZ], expectedZResult, 0.001_mm);
  }
}

}  // namespace Test
}  // namespace Acts
