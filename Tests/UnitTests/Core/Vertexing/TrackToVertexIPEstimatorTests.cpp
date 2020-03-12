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
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/TrackToVertexIPEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<> vXYDist(-0.1_mm, 0.1_mm);
// Vertex z position distribution
std::uniform_real_distribution<> vZDist(-20_mm, 20_mm);
// Track d0 distribution
std::uniform_real_distribution<> d0Dist(-0.01_mm, 0.01_mm);
// Track z0 distribution
std::uniform_real_distribution<> z0Dist(-0.2_mm, 0.2_mm);
// Track pT distribution
std::uniform_real_distribution<> pTDist(0.4_GeV, 10_GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100_um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Number of tracks distritbution
std::uniform_int_distribution<> nTracksDist(3, 10);

///
/// @brief Unit test for TrackToVertexIPEstimatorTests
///
BOOST_AUTO_TEST_CASE(track_to_vertex_ip_estimator_test) {
  // Number of tracks to test with
  unsigned int nTracks = 10;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  ConstantBField bField(0.0, 0.0, 1_T);

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator =
      std::make_shared<Propagator<EigenStepper<ConstantBField>>>(stepper);

  PropagatorOptions<> pOptions(tgContext, mfContext);

  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  // Create position of vertex and perigee surface
  double x = vXYDist(gen);
  double y = vXYDist(gen);
  double z = vZDist(gen);

  SpacePointVector vertexPosition(x, y, z, 0.);

  // Constraint for vertex fit
  Vertex<BoundParameters> myConstraint;
  // Some abitrary values
  SpacePointSymMatrix myCovMat = SpacePointSymMatrix::Zero();
  myCovMat(0, 0) = 30.;
  myCovMat(1, 1) = 30.;
  myCovMat(2, 2) = 30.;
  myCovMat(3, 3) = 30.;
  myConstraint.setFullCovariance(std::move(myCovMat));
  myConstraint.setFullPosition(vertexPosition);

  // Calculate d0 and z0 corresponding to vertex position
  double d0_v = sqrt(x * x + y * y);
  double z0_v = z;

  using IPEst_t =
      TrackToVertexIPEstimator<BoundParameters,
                               Propagator<EigenStepper<ConstantBField>>>;

  IPEst_t::Config ipEstCfg(propagator, pOptions);

  // Create TrackToVertexIPEstimator
  IPEst_t ipEst(ipEstCfg);

  // Construct random track emerging from vicinity of vertex position
  // Vector to store track objects used for vertex fit
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Construct random track parameters
    BoundVector paramVec;
    paramVec << d0_v + d0Dist(gen), z0_v + z0Dist(gen), phiDist(gen),
        thetaDist(gen), q / pTDist(gen), 0.;

    // Resolutions
    double resD0 = resIPDist(gen);
    double resZ0 = resIPDist(gen);
    double resPh = resAngDist(gen);
    double resTh = resAngDist(gen);
    double resQp = resQoPDist(gen);

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;
    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

    BoundParameters track =
        BoundParameters(tgContext, std::move(covMat), paramVec, perigeeSurface);

    // Check if IP are retrieved
    ImpactParametersAndSigma output =
        ipEst.estimate(track, myConstraint).value();
    BOOST_CHECK_NE(output.IPd0, 0.);
    BOOST_CHECK_NE(output.IPz0, 0.);
  }
}
}  // namespace Test
}  // namespace Acts
