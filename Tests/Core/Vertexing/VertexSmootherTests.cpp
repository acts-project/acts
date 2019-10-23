// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE VertexSmoother Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include <Acts/Vertexing/VertexSmoother.hpp>
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer_t = HelicalTrackLinearizer<ConstantBField, Propagator>;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<> vXYDist(-0.1 * units::_mm, 0.1 * units::_mm);
// Vertex z position distribution
std::uniform_real_distribution<> vZDist(-20 * units::_mm, 20 * units::_mm);
// Track d0 distribution
std::uniform_real_distribution<> d0Dist(-0.01 * units::_mm, 0.01 * units::_mm);
// Track z0 distribution
std::uniform_real_distribution<> z0Dist(-0.2 * units::_mm, 0.2 * units::_mm);
// Track pT distribution
std::uniform_real_distribution<> pTDist(0.4 * units::_GeV, 10. * units::_GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100. * units::_um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Number of tracks distritbution
std::uniform_int_distribution<> nTracksDist(3, 10);

/// @brief Unit test for VertexSmoother
///
BOOST_AUTO_TEST_CASE(sequential_vertex_smoother_test) {
  bool debugMode = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Number of tracks
  unsigned int nTracks = nTracksDist(gen);

  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up LinearizedTrackFactory, needed for linearizing the tracks
  PropagatorOptions<ActionList<>, AbortList<>> pOptions =
      Linearizer_t::getDefaultPropagatorOptions(tgContext, mfContext);

  // Linearizer for BoundParameters type test
  Linearizer_t::Config ltConfig(bField, propagator, pOptions);
  Linearizer_t linearizer(ltConfig);

  // Set up Billoir Vertex Fitter
  FullBilloirVertexFitter<BoundParameters, Linearizer_t>::Config
      vertexFitterCfg;
  FullBilloirVertexFitter<BoundParameters, Linearizer_t> billoirFitter(
      vertexFitterCfg);

  VertexFitterOptions<BoundParameters> vfOptions(tgContext, mfContext);

  // Now: create some tracks, fit and retrieve vertex using Billoirfitter
  // Create position of vertex and perigee surface
  double x = vXYDist(gen);
  double y = vXYDist(gen);
  double z = vZDist(gen);

  SpacePointVector vertexPosition(x, y, z, 0.);
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

  // Calculate d0 and z0 corresponding to vertex position
  double d0V = sqrt(x * x + y * y);
  double z0V = z;

  // Start constructing nTracks tracks in the following
  std::vector<BoundParameters> tracks;

  // Construct random track emerging from vicinity of vertex position
  // Vector to store track objects used for vertex fit
  for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Construct random track parameters
    BoundParameters::ParVector_t paramVec;
    paramVec << d0V + d0Dist(gen), z0V + z0Dist(gen), phiDist(gen),
        thetaDist(gen), q / pTDist(gen), 0.;

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;

    // Resolutions
    double resD0 = resIPDist(gen);
    double resZ0 = resIPDist(gen);
    double resPh = resAngDist(gen);
    double resTh = resAngDist(gen);
    double resQp = resQoPDist(gen);

    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;
    tracks.push_back(BoundParameters(tgContext, std::move(covMat), paramVec,
                                     perigeeSurface));
  }

  Vertex<BoundParameters> fittedVertex =
      billoirFitter.fit(tracks, linearizer, vfOptions).value();

  // copy vertex for later comparison
  Vertex<BoundParameters> vertexBeforeSmoothing = fittedVertex;

  VertexSmoothing::smoothVertexSequentially<BoundParameters>(tgContext,
                                                             &fittedVertex);

  // Billoirfitter does not provide the TracksAtVertex with a linearized
  // state,
  // hence returns tracks not refitted. However, sizes should match.
  BOOST_CHECK_EQUAL(vertexBeforeSmoothing.tracks().size(),
                    fittedVertex.tracks().size());

  std::vector<TrackAtVertex<BoundParameters>> tracksWithLinState;
  for (auto trackAtVtx : fittedVertex.tracks()) {
    BoundParameters fittedParams = trackAtVtx.fittedParams;

    LinearizedTrack linTrack =
        linearizer.linearizeTrack(&fittedParams, vertexPosition).value();
    trackAtVtx.linearizedState = linTrack;
    tracksWithLinState.push_back(trackAtVtx);
  }

  // set tracks with linearized state to vertex
  fittedVertex.setTracksAtVertex(tracksWithLinState);
  VertexSmoothing::smoothVertexSequentially<BoundParameters>(tgContext,
                                                             &fittedVertex);

  BOOST_CHECK_EQUAL(vertexBeforeSmoothing.tracks().size(),
                    fittedVertex.tracks().size());

  for (unsigned int i = 0; i < fittedVertex.tracks().size(); ++i) {
    auto paramOld = vertexBeforeSmoothing.tracks()[i].fittedParams;
    auto paramNew = fittedVertex.tracks()[i].fittedParams;
    BOOST_CHECK_NE(paramOld, paramNew);
    if (debugMode) {
      std::cout << "Track %d, old params: " << paramOld << std::endl;
      std::cout << "Track %d, new params: " << paramNew << std::endl;
    }
  }
}

}  // namespace Test
}  // namespace Acts
