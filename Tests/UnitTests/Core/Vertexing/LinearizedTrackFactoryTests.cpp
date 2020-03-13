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
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Linearizer =
    HelicalTrackLinearizer<Propagator<EigenStepper<ConstantBField>>>;

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

///
/// @brief Unit test for HelicalTrackLinearizer
///
BOOST_AUTO_TEST_CASE(linearized_track_factory_test) {
  // Number of tracks
  unsigned int nTracks = 100;

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

  // Calculate d0 and z0 corresponding to vertex position
  double d0v = sqrt(x * x + y * y);
  double z0v = z;

  // Start constructing nTracks tracks in the following
  std::vector<BoundParameters> tracks;

  // Construct random track emerging from vicinity of vertex position
  // Vector to store track objects used for vertex fit
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Construct random track parameters
    BoundVector paramVec;
    paramVec << d0v + d0Dist(gen), z0v + z0Dist(gen), phiDist(gen),
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
    tracks.push_back(BoundParameters(tgContext, std::move(covMat), paramVec,
                                     perigeeSurface));
  }

  Linearizer::Config ltConfig(bField, propagator, pOptions);
  Linearizer linFactory(ltConfig);

  BoundVector vecBoundZero = BoundVector::Zero();
  BoundSymMatrix matBoundZero = BoundSymMatrix::Zero();
  SpacePointVector vecSPZero = SpacePointVector::Zero();
  SpacePointToBoundMatrix matBound2SPZero = SpacePointToBoundMatrix::Zero();
  ActsMatrixD<BoundParsDim, 3> matBound2MomZero =
      ActsMatrixD<BoundParsDim, 3>::Zero();

  for (const BoundParameters& parameters : tracks) {
    LinearizedTrack linTrack =
        linFactory.linearizeTrack(parameters, SpacePointVector::Zero()).value();

    BOOST_CHECK_NE(linTrack.parametersAtPCA, vecBoundZero);
    BOOST_CHECK_NE(linTrack.covarianceAtPCA, matBoundZero);
    BOOST_CHECK_EQUAL(linTrack.linearizationPoint, vecSPZero);
    BOOST_CHECK_NE(linTrack.positionJacobian, matBound2SPZero);
    BOOST_CHECK_NE(linTrack.momentumJacobian, matBound2MomZero);
    BOOST_CHECK_NE(linTrack.constantTerm, vecBoundZero);
  }
}

///
/// @brief Unit test for HelicalTrackLinearizer
///
BOOST_AUTO_TEST_CASE(linearized_track_factory_straightline_test) {
  using LinearizerStraightLine =
      HelicalTrackLinearizer<Propagator<StraightLineStepper>>;
  // Number of tracks
  unsigned int nTracks = 100;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up stepper
  StraightLineStepper stepper;

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator<StraightLineStepper>>(stepper);

  PropagatorOptions<> pOptions(tgContext, mfContext);
  // Create perigee surface
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  // Create position of vertex and perigee surface
  double x = vXYDist(gen);
  double y = vXYDist(gen);
  double z = vZDist(gen);

  // Calculate d0 and z0 corresponding to vertex position
  double d0v = sqrt(x * x + y * y);
  double z0v = z;

  // Start constructing nTracks tracks in the following
  std::vector<BoundParameters> tracks;

  // Construct random track emerging from vicinity of vertex position
  // Vector to store track objects used for vertex fit
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Construct random track parameters
    BoundVector paramVec;
    paramVec << d0v + d0Dist(gen), z0v + z0Dist(gen), phiDist(gen),
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
    tracks.push_back(BoundParameters(tgContext, std::move(covMat), paramVec,
                                     perigeeSurface));
  }

  // Set up helical track linearizer for the case of a non-existing
  // magnetic field, which results in the extreme case of a straight line
  LinearizerStraightLine::Config ltConfig(propagator, pOptions);
  LinearizerStraightLine linFactory(ltConfig);

  BoundVector vecBoundZero = BoundVector::Zero();
  BoundSymMatrix matBoundZero = BoundSymMatrix::Zero();
  SpacePointVector vecSPZero = SpacePointVector::Zero();
  SpacePointToBoundMatrix matBound2SPZero = SpacePointToBoundMatrix::Zero();
  ActsMatrixD<BoundParsDim, 3> matBound2MomZero =
      ActsMatrixD<BoundParsDim, 3>::Zero();

  for (const BoundParameters& parameters : tracks) {
    LinearizedTrack linTrack =
        linFactory.linearizeTrack(parameters, SpacePointVector::Zero()).value();

    BOOST_CHECK_NE(linTrack.parametersAtPCA, vecBoundZero);
    BOOST_CHECK_NE(linTrack.covarianceAtPCA, matBoundZero);
    BOOST_CHECK_EQUAL(linTrack.linearizationPoint, vecSPZero);
    BOOST_CHECK_NE(linTrack.positionJacobian, matBound2SPZero);
    BOOST_CHECK_NE(linTrack.momentumJacobian, matBound2MomZero);
    BOOST_CHECK_NE(linTrack.constantTerm, vecBoundZero);
  }
}

}  // namespace Test
}  // namespace Acts
