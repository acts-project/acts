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

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

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
std::uniform_real_distribution<> pTDist(1._GeV, 30._GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100._um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.1, 0.1);
// Number of tracks distritbution
std::uniform_int_distribution<> nTracksDist(3, 10);

/// @brief Unit test for AdaptiveMultiVertexFitter
///
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_fitter_test) {
  bool debugMode = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 1._T));

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPoint3dEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator, pOptions);
  IPEstimator ip3dEst(ip3dEstCfg);

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer>::Config fitterCfg(
      ip3dEst);

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator, pOptions);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer> fitter(fitterCfg);

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer>::State state;

  // Create positions of three vertices, two of which (1 and 2) are
  // close to one another and will share a common track later
  Vector3D vtxPos1(-0.15_mm, -0.1_mm, -1.5_mm);
  Vector3D vtxPos2(-0.1_mm, -0.15_mm, -3._mm);
  Vector3D vtxPos3(0.2_mm, 0.2_mm, 10._mm);

  std::vector<Vector3D> vtxVec{vtxPos1, vtxPos2, vtxPos3};

  // Vector to store vectors of Tracks at vertex for every vertex
  std::vector<std::vector<TrackAtVertex<BoundParameters>>> trackVtxVec(
      vtxVec.size());

  // only for debugging
  std::vector<TrackAtVertex<BoundParameters>> allTracks;

  // Resolutions, use the same for all tracks
  double resD0 = resIPDist(gen);
  double resZ0 = resIPDist(gen);
  double resPh = resAngDist(gen);
  double resTh = resAngDist(gen);
  double resQp = resQoPDist(gen);

  unsigned int nTracksPerVtx = 4;
  // Construct nTracksPerVtx * 3 (3 vertices) random track emerging
  // from vicinity of vertex positions
  for (unsigned int iTrack = 0; iTrack < nTracksPerVtx * vtxVec.size();
       iTrack++) {
    // Construct positive or negative charge randomly
    double q = qDist(gen) < 0 ? -1. : 1.;

    // Fill vector of track objects with simple covariance matrix
    Covariance covMat;

    covMat << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0., 0.,
        0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh, 0.,
        0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

    // Index of current vertex
    int vtxIdx = (int)(iTrack / nTracksPerVtx);

    // Construct random track parameters
    BoundParameters::ParVector_t paramVec;
    paramVec << d0Dist(gen), z0Dist(gen), phiDist(gen), thetaDist(gen),
        q / pTDist(gen), 0.;

    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(vtxVec[vtxIdx]);

    auto trk =
        BoundParameters(tgContext, std::move(covMat), paramVec, perigeeSurface);

    TrackAtVertex<BoundParameters> trkAtVtx(1., trk, trk);

    if (debugMode) {
      std::cout << "Adding track " << iTrack << " for vertex " << vtxIdx
                << "\n\twith ID: " << trkAtVtx.id
                << "\n\tparams:  " << trk.parameters() << std::endl;
      allTracks.push_back(trkAtVtx);
    }

    trackVtxVec[vtxIdx].push_back(trkAtVtx);

    // Use first track also for second vertex to let vtx1 and vtx2
    // share this track
    if (iTrack == 0) {
      trackVtxVec[1].push_back(trkAtVtx);
    }
  }

  std::vector<Vertex<BoundParameters>> vtxList;

  int idx = 0;
  for (auto& vtxPos : vtxVec) {
    Vertex<BoundParameters> vtx(vtxPos);
    // Set track for current vertex
    vtx.setTracksAtVertex(trackVtxVec[idx]);
    // Set some vertex covariance
    SpacePointSymMatrix posCovariance(SpacePointSymMatrix::Identity());
    vtx.setFullCovariance(posCovariance);
    // Add to vertex list
    vtxList.push_back(vtx);
    idx++;
  }

  for (auto& vtx : vtxList) {
    std::cout << &vtx << std::endl;
  }

  for (auto& vtx : vtxList) {
    // Add vertex link to each track
    for (auto& trkAtVtx : vtx.tracks()) {
      state.trkInfoMap[trkAtVtx.id].linksToVertices.push_back(&vtx);
    }

    if (debugMode) {
      std::cout << "Number of tracks at vertex " << &vtx << ": "
                << vtx.tracks().size() << std::endl;
    }
  }

  if (debugMode) {
    for (auto& trkAtVtx : allTracks) {
      auto links = state.trkInfoMap[trkAtVtx.id].linksToVertices;
      for (auto vtxLink : links) {
        std::cout << "Track with ID: " << trkAtVtx.id << " used by vertex "
                  << vtxLink << std::endl;
      }
    }
  }

  // Copy vertex seeds from state.vertexCollection to new
  // list in order to be able to compare later
  std::vector<Vertex<BoundParameters>> seedListCopy = vtxList;

  auto res1 = fitter.fit(state, vtxList[0], linearizer, fitterOptions);

  BOOST_CHECK(res1.ok());

  if (debugMode) {
    std::cout << "Vertex positions after fit of vertex 1 and 2:" << std::endl;
    std::cout << "Vtx 1, seed position:\n " << seedListCopy[0].fullPosition()
              << "\nFitted position:\n " << vtxList[0].fullPosition()
              << std::endl;
    std::cout << "Vtx 2, seed position:\n " << seedListCopy[1].fullPosition()
              << "\nFitted position:\n " << vtxList[1].fullPosition()
              << std::endl;
    std::cout << "Vtx 3, seed position:\n " << seedListCopy[2].fullPosition()
              << "\nFitted position:\n " << vtxList[2].fullPosition()
              << std::endl;
  }

  // After fit of first vertex, only first and second vertex seed
  // should have been modified while third vertex should remain untouched
  BOOST_CHECK_NE(vtxList[0].fullPosition(), seedListCopy[0].fullPosition());
  BOOST_CHECK_NE(vtxList[1].fullPosition(), seedListCopy[1].fullPosition());
  BOOST_CHECK_EQUAL(vtxList[2].fullPosition(), seedListCopy[2].fullPosition());

  CHECK_CLOSE_ABS(vtxList[0].fullPosition(), seedListCopy[0].fullPosition(),
                  1_mm);
  CHECK_CLOSE_ABS(vtxList[1].fullPosition(), seedListCopy[1].fullPosition(),
                  1_mm);

  auto res2 = fitter.fit(state, vtxList[2], linearizer, fitterOptions);

  BOOST_CHECK(res2.ok());

  // Now also the third vertex should have been modified and fitted
  BOOST_CHECK_NE(vtxList[2].fullPosition(), seedListCopy[2].fullPosition());
  CHECK_CLOSE_ABS(vtxList[2].fullPosition(), seedListCopy[2].fullPosition(),
                  1_mm);

  if (debugMode) {
    std::cout << "Vertex positions after fit of vertex 3:" << std::endl;
    std::cout << "Vtx 1, seed position:\n " << seedListCopy[0].fullPosition()
              << "\nFitted position:\n " << vtxList[0].fullPosition()
              << std::endl;
    std::cout << "Vtx 2, seed position:\n " << seedListCopy[1].fullPosition()
              << "\nFitted position:\n " << vtxList[1].fullPosition()
              << std::endl;
    std::cout << "Vtx 3, seed position:\n " << seedListCopy[2].fullPosition()
              << "\nFitted position:\n " << vtxList[2].fullPosition()
              << std::endl;
  }
}

/// @brief Unit test for AdaptiveMultiVertexFitter
/// based on Athena unit test, i.e. same setting and
/// test values are used here
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_fitter_test_athena) {
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2._T));

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPoint3dEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator, pOptions);
  IPEstimator ip3dEst(ip3dEstCfg);

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer>::Config fitterCfg(
      ip3dEst);

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator, pOptions);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  // fitterCfg.doSmoothing = true;

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer> fitter(fitterCfg);

  AdaptiveMultiVertexFitter<BoundParameters, Linearizer>::State state;

  // Create first vector of tracks
  Vector3D pos0(0., 0., 0.);
  Vector3D pos1a(2_mm, 1_mm, -10_mm);
  Vector3D mom1a(400_MeV, 600_MeV, 200_MeV);
  Vector3D pos1b(1_mm, 2_mm, -3_mm);
  Vector3D mom1b(600_MeV, 400_MeV, -200_MeV);
  Vector3D pos1c(1.2_mm, 1.3_mm, -7_mm);
  Vector3D mom1c(300_MeV, 1000_MeV, 100_MeV);

  // Start creating some track parameters
  Covariance covMat1 = Covariance::Identity();

  std::vector<BoundParameters> params1;

  params1.push_back(
      BoundParameters(tgContext, covMat1, pos1a, mom1a, 1, 0,
                      Surface::makeShared<PerigeeSurface>(pos1a)));
  params1.push_back(
      BoundParameters(tgContext, covMat1, pos1b, mom1b, -1, 0,
                      Surface::makeShared<PerigeeSurface>(pos1b)));
  params1.push_back(
      BoundParameters(tgContext, covMat1, pos1c, mom1c, -1, 0,
                      Surface::makeShared<PerigeeSurface>(pos1c)));

  // Create second vector of tracks
  Vector3D pos2a(10_mm, 0_mm, -5_mm);
  Vector3D mom2a(1000_MeV, 0_MeV, 0_MeV);
  Vector3D pos2b(10.5_mm, 0.5_mm, -5.5_mm);
  Vector3D mom2b(800_MeV, 200_MeV, 200_MeV);
  Vector3D pos2c(9.5_mm, -0.5_mm, -4.5_mm);
  Vector3D mom2c(700_MeV, -300_MeV, -200_MeV);
  Covariance covMat2 = Covariance::Identity();

  // Define covariance entries as used in athena unit test
  const double covEntries = 1e-2;
  covMat2 = covMat2 * covEntries;
  covMat2(1, 1) = 1;

  std::vector<BoundParameters> params2;

  params2.push_back(
      BoundParameters(tgContext, covMat2, pos2a, mom2a, 1, 0,
                      Surface::makeShared<PerigeeSurface>(pos2a)));
  params2.push_back(
      BoundParameters(tgContext, covMat2, pos2b, mom2b, -1, 0,
                      Surface::makeShared<PerigeeSurface>(pos2b)));
  params2.push_back(
      BoundParameters(tgContext, covMat2, pos2c, mom2c, -1, 0,
                      Surface::makeShared<PerigeeSurface>(pos2c)));

  std::vector<TrackAtVertex<BoundParameters>> tracksAtVtx1;
  for (const auto& trk : params1) {
    tracksAtVtx1.push_back(TrackAtVertex<BoundParameters>(0, trk, trk));
  }
  std::vector<TrackAtVertex<BoundParameters>> tracksAtVtx2;
  for (const auto& trk : params2) {
    tracksAtVtx2.push_back(TrackAtVertex<BoundParameters>(0, trk, trk));
  }

  std::vector<Vertex<BoundParameters>> vtxList;

  Vector3D vtxPos1(1.5_mm, 1.7_mm, -6_mm);
  Vertex<BoundParameters> vtx1(vtxPo1);
  // Set track for current vertex
  vtx1.setTracksAtVertex(tracksAtVtx1);
  // Add to vertex list
  vtxList.push_back(vtx1);

  Vector3D vtxPos2(9.8_mm, 0.2_mm, -4.8_mm);
  Vertex<BoundParameters> vtx2(vtxPo1);
  // Set track for current vertex
  vtx2.setTracksAtVertex(tracksAtVtx2);
  // Add to vertex list
  vtxList.push_back(vtx2);
}

}  // namespace Test
}  // namespace Acts
