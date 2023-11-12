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
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/KalmanVertexTrackUpdater.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSquareMatrix;
using Propagator = Acts::Propagator<EigenStepper<>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

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
std::uniform_real_distribution<> resQoPDist(-0.01, 0.01);

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
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;
  IPEstimator::Config ip3dEstConfig(bField, propagator);
  IPEstimator ip3dEst(ip3dEstConfig);
  IPEstimator::State state(bField->makeCache(magFieldContext));

  // Set up HelicalTrackLinearizer, needed for linearizing the tracks
  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);
  Linearizer::State linState(bField->makeCache(magFieldContext));

  // Create perigee surface at origin
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3(0., 0., 0.));

  // Create random tracks around origin and a random vertex.
  // Update tracks with the assumption that they originate from
  // the vertex position and check if they are closer to the
  // vertex after the update process
  for (unsigned int i = 0; i < nTests; ++i) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Construct random track parameters
    BoundTrackParameters::ParametersVector paramVec;

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
        linearizer
            .linearizeTrack(params, 0, *perigee, geoContext, magFieldContext,
                            linState)
            .value();

    // Create TrackAtVertex
    TrackAtVertex<BoundTrackParameters> trkAtVtx(0., params, &params);

    // Set linearized state of trackAtVertex
    trkAtVtx.linearizedState = linTrack;

    // Copy parameters for later comparison of old and new version
    auto fittedParamsCopy = trkAtVtx.fittedParams;

    // Create a vertex
    Vector3 vtxPos(vXYDist(gen), vXYDist(gen), vZDist(gen));
    Vertex<BoundTrackParameters> vtx(vtxPos);

    // Update trkAtVertex with assumption of originating from vtx
    KalmanVertexTrackUpdater::update<BoundTrackParameters>(trkAtVtx, vtx);

    // The old distance
    double oldDistance =
        ip3dEst.calculateDistance(geoContext, fittedParamsCopy, vtxPos, state)
            .value();

    // The new distance after update
    double newDistance =
        ip3dEst
            .calculateDistance(geoContext, trkAtVtx.fittedParams, vtxPos, state)
            .value();
    if (debug) {
      std::cout << "Old distance: " << oldDistance << std::endl;
      std::cout << "New distance: " << newDistance << std::endl;
    }

    // Parameters should have changed
    BOOST_CHECK_NE(fittedParamsCopy, trkAtVtx.fittedParams);

    // After update, track should be closer to the vertex
    BOOST_CHECK_LT(newDistance, oldDistance);

  }  // end for loop

}  // end test case

}  // namespace Test
}  // namespace Acts
