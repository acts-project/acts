// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldQuery.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/detail/KalmanVertexUpdaterImpl.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using Covariance = BoundMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer;

// Create a test context
GeometryContext geoContext = GeometryContext::dangerouslyDefaultConstruct();
MagneticFieldContext magFieldContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<double> vXYDist(-0.1_mm, 0.1_mm);
// Vertex z position distribution
std::uniform_real_distribution<double> vZDist(-20_mm, 20_mm);
// Track d0 distribution
std::uniform_real_distribution<double> d0Dist(-0.01_mm, 0.01_mm);
// Track z0 distribution
std::uniform_real_distribution<double> z0Dist(-0.2_mm, 0.2_mm);
// Track pT distribution
std::uniform_real_distribution<double> pTDist(0.4_GeV, 10_GeV);
// Track phi distribution
std::uniform_real_distribution<double> phiDist(-std::numbers::pi,
                                               std::numbers::pi);
// Track theta distribution
std::uniform_real_distribution<double> thetaDist(1., std::numbers::pi - 1.);
// Track charge helper distribution
std::uniform_real_distribution<double> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<double> resIPDist(0., 100_um);
// Track angular distribution
std::uniform_real_distribution<double> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<double> resQoPDist(-0.01, 0.01);
// Number of vertices per test event distribution

BOOST_AUTO_TEST_SUITE(VertexingSuite)
///
/// @brief Unit test for KalmanVertexUpdater
///
BOOST_AUTO_TEST_CASE(Kalman_Vertex_Updater) {
  bool debug = false;

  // Number of tests
  unsigned int nTests = 10;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig;

  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);
  MagneticFieldQuery fieldQuery(*bField, magFieldContext);

  // Create perigee surface at origin
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  // Creates a random tracks around origin and a random vertex.
  // VertexUpdater adds track to vertex and updates the position
  // which should afterwards be closer to the origin/track
  for (unsigned int i = 0; i < nTests; ++i) {
    if (debug) {
      std::cout << "Test " << i + 1 << std::endl;
    }
    // Construct positive or negative charge randomly
    double q = std::copysign(1., qDist(gen));

    // Construct random track parameters around origin
    BoundVector paramVec;

    paramVec << d0Dist(gen), z0Dist(gen), phiDist(gen), thetaDist(gen),
        q / pTDist(gen), 0.;

    if (debug) {
      std::cout << "Creating track parameters: " << paramVec << std::endl;
    }

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;

    // Resolutions
    double res_d0 = resIPDist(gen);
    double res_z0 = resIPDist(gen);
    double res_ph = resAngDist(gen);
    double res_th = resAngDist(gen);
    double res_qp = resQoPDist(gen);

    covMat << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0, 0., 0.,
        0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
        res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0., 0.,
        0., 0., 0., 1.;
    BoundTrackParameters params(perigeeSurface, paramVec, std::move(covMat),
                                ParticleHypothesis::pion());

    std::shared_ptr<PerigeeSurface> perigee =
        Surface::makeShared<PerigeeSurface>(Vector3::Zero());

    // Linearized state of the track
    LinearizedTrack linTrack =
        linearizer.linearizeTrack(params, 0, *perigee, geoContext, fieldQuery)
            .value();

    // Create TrackAtVertex
    TrackAtVertex trkAtVtx(0., params, InputTrack{&params});

    // Set linearized state of trackAtVertex
    trkAtVtx.linearizedState = linTrack;
    trkAtVtx.isLinearized = true;

    // Create a vertex position
    Vector3 vtxPos(vXYDist(gen), vXYDist(gen), vZDist(gen));

    // Fit the spatial vertex coordinates (3) or also time (4)
    for (unsigned int nDimVertex : {3u, 4u}) {
      Vertex vtx(vtxPos);
      vtx.setFullCovariance(SquareMatrix4::Identity() * 0.01);

      // Update trkAtVertex with assumption of originating from vtx
      KalmanVertexUpdater::updateVertexWithTrack(vtx, trkAtVtx, nDimVertex);

      if (debug) {
        std::cout << "Old vertex position: " << vtxPos << std::endl;
        std::cout << "New vertex position: " << vtx.position() << std::endl;
      }

      double oldDistance = vtxPos.norm();
      double newDistance = vtx.position().norm();

      if (debug) {
        std::cout << "Old distance: " << oldDistance << std::endl;
        std::cout << "New distance: " << newDistance << std::endl;
      }

      // After update, vertex should be closer to the track
      BOOST_CHECK_LT(newDistance, oldDistance);

      // Note: KalmanVertexUpdater updates the vertex w.r.t. the
      // newly given track, but does NOT add the track to the
      // TrackAtVertex list. Has to be done manually after calling
      // the update method.
      BOOST_CHECK(vtx.tracks().empty());
    }

  }  // end for loop

}  // end test case

///
/// @brief Unit test for KalmanVertexTrackUpdater
///
BOOST_AUTO_TEST_CASE(Kalman_Vertex_TrackUpdater) {
  bool debug = true;

  // Number of tests
  unsigned int nTests = 10;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up Eigenstepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up ImpactPointEstimator, used for comparisons later
  ImpactPointEstimator::Config ip3dEstConfig(bField, propagator);
  ImpactPointEstimator ip3dEst(ip3dEstConfig);
  ImpactPointEstimator::State state{bField->makeCache(magFieldContext)};

  // Set up HelicalTrackLinearizer, needed for linearizing the tracks
  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig;

  ltConfig.propagator = propagator;
  Linearizer linearizer(ltConfig);
  MagneticFieldQuery fieldQuery(*bField, magFieldContext);

  // Create perigee surface at origin
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  // Create random tracks around origin and a random vertex.
  // Update tracks with the assumption that they originate from
  // the vertex position and check if they are closer to the
  // vertex after the update process
  for (unsigned int i = 0; i < nTests; ++i) {
    // Construct positive or negative charge randomly
    double q = std::copysign(1., qDist(gen));

    // Construct random track parameters
    BoundVector paramVec;

    paramVec << d0Dist(gen), z0Dist(gen), phiDist(gen), thetaDist(gen),
        q / pTDist(gen), 0.;

    if (debug) {
      std::cout << "Creating track parameters: " << paramVec << std::endl;
    }

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;

    // Resolutions
    double res_d0 = resIPDist(gen);
    double res_z0 = resIPDist(gen);
    double res_ph = resAngDist(gen);
    double res_th = resAngDist(gen);
    double res_qp = resQoPDist(gen);

    covMat << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0, 0., 0.,
        0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
        res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0., 0.,
        0., 0., 0., 1.;
    BoundTrackParameters params(perigeeSurface, paramVec, std::move(covMat),
                                ParticleHypothesis::pion());

    std::shared_ptr<PerigeeSurface> perigee =
        Surface::makeShared<PerigeeSurface>(Vector3::Zero());

    // Linearized state of the track
    LinearizedTrack linTrack =
        linearizer.linearizeTrack(params, 0, *perigee, geoContext, fieldQuery)
            .value();

    // Create a vertex position
    Vector3 vtxPos(vXYDist(gen), vXYDist(gen), vZDist(gen));

    // Fit the spatial vertex coordinates (3) or also time (4)
    for (unsigned int nDimVertex : {3u, 4u}) {
      // Create TrackAtVertex
      TrackAtVertex trkAtVtx(0., params, InputTrack{&params});

      // Set linearized state of trackAtVertex
      trkAtVtx.linearizedState = linTrack;
      trkAtVtx.isLinearized = true;

      // Copy parameters for later comparison of old and new version
      auto fittedParamsCopy = trkAtVtx.fittedParams;

      Vertex vtx(vtxPos);

      // Update trkAtVertex with assumption of originating from vtx
      KalmanVertexUpdater::updateTrackWithVertex(trkAtVtx, vtx, nDimVertex);

      // The old distance
      double oldDistance =
          ip3dEst.calculateDistance(geoContext, fittedParamsCopy, vtxPos, state)
              .value();

      // The new distance after update
      double newDistance =
          ip3dEst
              .calculateDistance(geoContext, trkAtVtx.fittedParams, vtxPos,
                                 state)
              .value();
      if (debug) {
        std::cout << "Old distance: " << oldDistance << std::endl;
        std::cout << "New distance: " << newDistance << std::endl;
      }

      // Parameters should have changed
      BOOST_CHECK_NE(fittedParamsCopy, trkAtVtx.fittedParams);

      // After update, track should be closer to the vertex
      BOOST_CHECK_LT(newDistance, oldDistance);
    }

  }  // end for loop

}  // end test case

///
/// @brief Unit test for the weight of the fitted track parameters: the full
/// weight when fitting time, otherwise the inverse of the covariance without
/// time, obtained from the full weight or, without time information, from the
/// covariance
///
BOOST_AUTO_TEST_CASE(Kalman_Vertex_TrackParameterWeight) {
  using WeightWithoutTime = SquareMatrix<eBoundSize - 1>;
  std::mt19937 gen(31415);
  std::normal_distribution<double> gauss(0., 1.);

  // Random covariance with the typical scales of perigee parameters (mm, rad,
  // 1/GeV, ns), with the time correlations scaled by `timeCorrelationScale`
  auto makeCovariance = [&](double timeSigma, double timeCorrelationScale) {
    BoundMatrix a;
    for (unsigned int i = 0; i < eBoundSize; ++i) {
      for (unsigned int j = 0; j < eBoundSize; ++j) {
        a(i, j) = gauss(gen);
      }
    }
    const BoundMatrix spd = a * a.transpose() + 6. * BoundMatrix::Identity();
    const BoundVector invSigma = spd.diagonal().cwiseSqrt().cwiseInverse();
    BoundMatrix correlation =
        invSigma.asDiagonal() * spd * invSigma.asDiagonal();
    correlation.row(eBoundTime).head<eBoundSize - 1>() *= timeCorrelationScale;
    correlation.col(eBoundTime).head<eBoundSize - 1>() *= timeCorrelationScale;
    BoundVector sigma;
    sigma << 20_um, 50_um, 1e-4, 1e-4, 1e-3 / 1_GeV, timeSigma;
    return BoundMatrix(sigma.asDiagonal() * correlation * sigma.asDiagonal());
  };

  struct TimeCase {
    double sigma;
    double correlationScale;
  };
  // Correlated time, uncorrelated time, precise and inflated time
  for (const TimeCase& timeCase :
       {TimeCase{1_ns, 1.}, TimeCase{1_ns, 0.}, TimeCase{1e-3_ns, 1.},
        TimeCase{1e6_ns, 1.}}) {
    for (unsigned int i = 0; i < 100; ++i) {
      LinearizedTrack linTrack;
      linTrack.covarianceAtPCA =
          makeCovariance(timeCase.sigma, timeCase.correlationScale);
      linTrack.weightAtPCA = linTrack.covarianceAtPCA.inverse();

      // Fitting time: the full weight
      BOOST_CHECK_EQUAL(
          KalmanVertexUpdater::detail::trackParameterWeight<eBoundSize>(
              linTrack),
          linTrack.weightAtPCA);

      // Not fitting time: the inverse of the covariance without time
      const WeightWithoutTime weight =
          KalmanVertexUpdater::detail::trackParameterWeight<eBoundSize - 1>(
              linTrack);
      const WeightWithoutTime covarianceWithoutTime =
          linTrack.covarianceAtPCA
              .topLeftCorner<eBoundSize - 1, eBoundSize - 1>();
      CHECK_CLOSE_OR_SMALL(weight, covarianceWithoutTime.inverse().eval(),
                           1e-10, 1e-10);
      CHECK_CLOSE_OR_SMALL((weight * covarianceWithoutTime).eval(),
                           WeightWithoutTime::Identity(), 1e-10, 1e-10);
    }
  }

  // Without time information the full covariance is singular, and the
  // covariance without time is inverted directly
  for (unsigned int i = 0; i < 100; ++i) {
    LinearizedTrack linTrack;
    linTrack.covarianceAtPCA = makeCovariance(0., 0.);
    linTrack.weightAtPCA = linTrack.covarianceAtPCA.inverse();
    const WeightWithoutTime covarianceWithoutTime =
        linTrack.covarianceAtPCA
            .topLeftCorner<eBoundSize - 1, eBoundSize - 1>();
    const WeightWithoutTime weight =
        KalmanVertexUpdater::detail::trackParameterWeight<eBoundSize - 1>(
            linTrack);
    BOOST_CHECK(weight.allFinite());
    BOOST_CHECK_EQUAL(weight, covarianceWithoutTime.inverse().eval());
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
