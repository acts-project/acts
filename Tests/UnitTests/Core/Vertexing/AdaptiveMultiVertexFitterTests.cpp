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
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

using Covariance = BoundSymMatrix;
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
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 1_T});

  // Set up EigenStepper
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator);
  IPEstimator ip3dEst(ip3dEstCfg);

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>::Config fitterCfg(
      ip3dEst);

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer> fitter(fitterCfg);

  // Create positions of three vertices, two of which (1 and 2) are
  // close to one another and will share a common track later
  Vector3 vtxPos1(-0.15_mm, -0.1_mm, -1.5_mm);
  Vector3 vtxPos2(-0.1_mm, -0.15_mm, -3._mm);
  Vector3 vtxPos3(0.2_mm, 0.2_mm, 10._mm);

  std::vector<Vector3> vtxPosVec{vtxPos1, vtxPos2, vtxPos3};

  // Resolutions, use the same for all tracks
  double resD0 = resIPDist(gen);
  double resZ0 = resIPDist(gen);
  double resPh = resAngDist(gen);
  double resTh = resAngDist(gen);
  double resQp = resQoPDist(gen);

  std::vector<Vertex<BoundTrackParameters>> vtxList;
  for (auto& vtxPos : vtxPosVec) {
    Vertex<BoundTrackParameters> vtx(vtxPos);
    // Set some vertex covariance
    SymMatrix4 posCovariance(SymMatrix4::Identity());
    vtx.setFullCovariance(posCovariance);
    // Add to vertex list
    vtxList.push_back(vtx);
  }

  std::vector<Vertex<BoundTrackParameters>*> vtxPtrList;
  int cv = 0;
  if (debugMode) {
    std::cout << "All vertices in test case: " << std::endl;
  }
  for (auto& vtx : vtxList) {
    if (debugMode) {
      cv++;
      std::cout << "\t" << cv << ". vertex ptr: " << &vtx << std::endl;
    }
    vtxPtrList.push_back(&vtx);
  }

  std::vector<BoundTrackParameters> allTracks;

  unsigned int nTracksPerVtx = 4;
  // Construct nTracksPerVtx * 3 (3 vertices) random track emerging
  // from vicinity of vertex positions
  for (unsigned int iTrack = 0; iTrack < nTracksPerVtx * vtxPosVec.size();
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
    BoundTrackParameters::ParametersVector paramVec;
    paramVec << d0Dist(gen), z0Dist(gen), phiDist(gen), thetaDist(gen),
        q / pTDist(gen), 0.;

    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(vtxPosVec[vtxIdx]);

    allTracks.emplace_back(perigeeSurface, paramVec, std::move(covMat));
  }

  if (debugMode) {
    int ct = 0;
    std::cout << "All tracks in test case: " << std::endl;
    for (auto& trk : allTracks) {
      ct++;
      std::cout << "\t" << ct << ". track ptr: " << &trk << std::endl;
    }
  }

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>::State state(
      *bField, magFieldContext);

  for (unsigned int iTrack = 0; iTrack < nTracksPerVtx * vtxPosVec.size();
       iTrack++) {
    // Index of current vertex
    int vtxIdx = (int)(iTrack / nTracksPerVtx);
    state.vtxInfoMap[&(vtxList[vtxIdx])].trackLinks.push_back(
        &(allTracks[iTrack]));
    state.tracksAtVerticesMap.insert(
        std::make_pair(std::make_pair(&(allTracks[iTrack]), &(vtxList[vtxIdx])),
                       TrackAtVertex<BoundTrackParameters>(
                           1., allTracks[iTrack], &(allTracks[iTrack]))));

    // Use first track also for second vertex to let vtx1 and vtx2
    // share this track
    if (iTrack == 0) {
      state.vtxInfoMap[&(vtxList.at(1))].trackLinks.push_back(
          &(allTracks[iTrack]));
      state.tracksAtVerticesMap.insert(
          std::make_pair(std::make_pair(&(allTracks[iTrack]), &(vtxList.at(1))),
                         TrackAtVertex<BoundTrackParameters>(
                             1., allTracks[iTrack], &(allTracks[iTrack]))));
    }
  }

  for (auto& vtx : vtxPtrList) {
    state.addVertexToMultiMap(*vtx);
    if (debugMode) {
      std::cout << "Vertex, with ptr: " << vtx << std::endl;
      for (auto& trk : state.vtxInfoMap[vtx].trackLinks) {
        std::cout << "\t track ptr: " << trk << std::endl;
      }
    }
  }

  if (debugMode) {
    std::cout << "Checking all vertices linked to a single track: "
              << std::endl;
    for (auto& trk : allTracks) {
      std::cout << "Track with ptr: " << &trk << std::endl;
      auto range = state.trackToVerticesMultiMap.equal_range(&trk);
      for (auto vtxIter = range.first; vtxIter != range.second; ++vtxIter) {
        std::cout << "\t used by vertex: " << vtxIter->second << std::endl;
      }
    }
  }

  // Copy vertex seeds from state.vertexCollection to new
  // list in order to be able to compare later
  std::vector<Vertex<BoundTrackParameters>> seedListCopy = vtxList;

  auto res1 =
      fitter.addVtxToFit(state, vtxList.at(0), linearizer, vertexingOptions);
  if (debugMode) {
    std::cout << "Tracks linked to each vertex AFTER fit: " << std::endl;
    int c = 0;
    for (auto& vtx : vtxPtrList) {
      c++;
      std::cout << c << ". vertex, with ptr: " << vtx << std::endl;
      for (auto& trk : state.vtxInfoMap[vtx].trackLinks) {
        std::cout << "\t track ptr: " << trk << std::endl;
      }
    }
  }

  if (debugMode) {
    std::cout << "Checking all vertices linked to a single track AFTER fit: "
              << std::endl;
    for (auto& trk : allTracks) {
      std::cout << "Track with ptr: " << &trk << std::endl;
      auto range = state.trackToVerticesMultiMap.equal_range(&trk);
      for (auto vtxIter = range.first; vtxIter != range.second; ++vtxIter) {
        std::cout << "\t used by vertex: " << vtxIter->second << std::endl;
      }
    }
  }

  BOOST_CHECK(res1.ok());

  if (debugMode) {
    std::cout << "Vertex positions after fit of vertex 1 and 2:" << std::endl;
    std::cout << "Vtx 1, seed position:\n " << seedListCopy.at(0).fullPosition()
              << "\nFitted position:\n " << vtxList.at(0).fullPosition()
              << std::endl;
    std::cout << "Vtx 2, seed position:\n " << seedListCopy.at(1).fullPosition()
              << "\nFitted position:\n " << vtxList.at(1).fullPosition()
              << std::endl;
    std::cout << "Vtx 3, seed position:\n " << seedListCopy.at(2).fullPosition()
              << "\nFitted position:\n " << vtxList.at(2).fullPosition()
              << std::endl;
  }

  // After fit of first vertex, only first and second vertex seed
  // should have been modified while third vertex should remain untouched
  BOOST_CHECK_NE(vtxList.at(0).fullPosition(),
                 seedListCopy.at(0).fullPosition());
  BOOST_CHECK_NE(vtxList.at(1).fullPosition(),
                 seedListCopy.at(1).fullPosition());
  BOOST_CHECK_EQUAL(vtxList.at(2).fullPosition(),
                    seedListCopy.at(2).fullPosition());

  CHECK_CLOSE_ABS(vtxList.at(0).fullPosition(),
                  seedListCopy.at(0).fullPosition(), 1_mm);
  CHECK_CLOSE_ABS(vtxList.at(1).fullPosition(),
                  seedListCopy.at(1).fullPosition(), 1_mm);

  auto res2 =
      fitter.addVtxToFit(state, vtxList.at(2), linearizer, vertexingOptions);
  BOOST_CHECK(res2.ok());

  // Now also the third vertex should have been modified and fitted
  BOOST_CHECK_NE(vtxList.at(2).fullPosition(),
                 seedListCopy.at(2).fullPosition());
  CHECK_CLOSE_ABS(vtxList.at(2).fullPosition(),
                  seedListCopy.at(2).fullPosition(), 1_mm);

  if (debugMode) {
    std::cout << "Vertex positions after fit of vertex 3:" << std::endl;
    std::cout << "Vtx 1, seed position:\n " << seedListCopy.at(0).fullPosition()
              << "\nFitted position:\n " << vtxList.at(0).fullPosition()
              << std::endl;
    std::cout << "Vtx 2, seed position:\n " << seedListCopy.at(1).fullPosition()
              << "\nFitted position:\n " << vtxList.at(1).fullPosition()
              << std::endl;
    std::cout << "Vtx 3, seed position:\n " << seedListCopy.at(2).fullPosition()
              << "\nFitted position:\n " << vtxList.at(2).fullPosition()
              << std::endl;
  }
}

/// @brief Unit test for AdaptiveMultiVertexFitter
/// based on Athena unit test, i.e. same setting and
/// test values are used here
BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_fitter_test_athena) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  auto bField = std::make_shared<ConstantBField>(Vector3{0.0, 0.0, 2_T});

  // Set up EigenStepper
  // EigenStepper<> stepper(bField);
  EigenStepper<> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                          magFieldContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundTrackParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator);
  IPEstimator ip3dEst(ip3dEstCfg);

  std::vector<double> temperatures(1, 3.);
  AnnealingUtility::Config annealingConfig;
  annealingConfig.setOfTemperatures = temperatures;
  AnnealingUtility annealingUtility(annealingConfig);

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>::Config fitterCfg(
      ip3dEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundTrackParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  // fitterCfg.doSmoothing = true;

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer> fitter(fitterCfg);

  // Create first vector of tracks
  Vector3 pos1a(0.5_mm, -0.5_mm, 2.4_mm);
  Vector3 mom1a(1000_MeV, 0_MeV, -500_MeV);
  Vector3 pos1b(0.5_mm, -0.5_mm, 3.5_mm);
  Vector3 mom1b(0_MeV, 1000_MeV, 500_MeV);
  Vector3 pos1c(-0.2_mm, 0.1_mm, 3.4_mm);
  Vector3 mom1c(-50_MeV, 180_MeV, 300_MeV);

  Vector3 pos1d(-0.1_mm, 0.3_mm, 3.0_mm);
  Vector3 mom1d(-80_MeV, 480_MeV, -100_MeV);
  Vector3 pos1e(-0.01_mm, 0.01_mm, 2.9_mm);
  Vector3 mom1e(-600_MeV, 10_MeV, 210_MeV);

  Vector3 pos1f(-0.07_mm, 0.03_mm, 2.5_mm);
  Vector3 mom1f(240_MeV, 110_MeV, 150_MeV);

  // Start creating some track parameters
  Covariance covMat1;
  covMat1 << 1_mm * 1_mm, 0, 0., 0, 0., 0, 0, 1_mm * 1_mm, 0, 0., 0, 0, 0., 0,
      0.1, 0, 0, 0, 0, 0., 0, 0.1, 0, 0, 0., 0, 0, 0, 1. / (10_GeV * 10_GeV), 0,
      0, 0, 0, 0, 0, 1_ns;

  std::vector<BoundTrackParameters> params1 = {
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1a),
                                   geoContext, makeVector4(pos1a, 0), mom1a,
                                   mom1a.norm(), 1, covMat1)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1b),
                                   geoContext, makeVector4(pos1b, 0), mom1b,
                                   mom1b.norm(), -1, covMat1)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1c),
                                   geoContext, makeVector4(pos1c, 0), mom1c,
                                   mom1c.norm(), 1, covMat1)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1d),
                                   geoContext, makeVector4(pos1d, 0), mom1d,
                                   mom1d.norm(), -1, covMat1)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1e),
                                   geoContext, makeVector4(pos1e, 0), mom1e,
                                   mom1e.norm(), 1, covMat1)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos1f),
                                   geoContext, makeVector4(pos1f, 0), mom1f,
                                   mom1f.norm(), -1, covMat1)
          .value(),
  };

  // Create second vector of tracks
  Vector3 pos2a(0.2_mm, 0_mm, -4.9_mm);
  Vector3 mom2a(5000_MeV, 30_MeV, 200_MeV);
  Vector3 pos2b(-0.5_mm, 0.1_mm, -5.1_mm);
  Vector3 mom2b(800_MeV, 1200_MeV, 200_MeV);
  Vector3 pos2c(0.05_mm, -0.5_mm, -4.7_mm);
  Vector3 mom2c(400_MeV, -300_MeV, -200_MeV);

  // Define covariance as used in athena unit test
  Covariance covMat2 = covMat1;

  std::vector<BoundTrackParameters> params2 = {
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos2a),
                                   geoContext, makeVector4(pos2a, 0), mom2a,
                                   mom2a.norm(), 1, covMat2)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos2b),
                                   geoContext, makeVector4(pos2b, 0), mom2b,
                                   mom2b.norm(), -1, covMat2)
          .value(),
      BoundTrackParameters::create(Surface::makeShared<PerigeeSurface>(pos2c),
                                   geoContext, makeVector4(pos2c, 0), mom2c,
                                   mom2c.norm(), -1, covMat2)
          .value(),
  };
  std::vector<Vertex<BoundTrackParameters>*> vtxList;

  AdaptiveMultiVertexFitter<BoundTrackParameters, Linearizer>::State state(
      *bField, magFieldContext);

  // The constraint vertex position covariance
  SymMatrix4 covConstr(SymMatrix4::Identity());
  covConstr = covConstr * 1e+8;
  covConstr(3, 3) = 0.;

  // Prepare first vertex
  Vector3 vtxPos1(0.15_mm, 0.15_mm, 2.9_mm);
  Vertex<BoundTrackParameters> vtx1(vtxPos1);

  // Add to vertex list
  vtxList.push_back(&vtx1);

  // The constraint vtx for vtx1
  Vertex<BoundTrackParameters> vtx1Constr(vtxPos1);
  vtx1Constr.setFullCovariance(covConstr);
  vtx1Constr.setFitQuality(0, -3);

  // Prepare vtx info for fitter
  VertexInfo<BoundTrackParameters> vtxInfo1;
  vtxInfo1.linPoint.setZero();
  vtxInfo1.linPoint.head<3>() = vtxPos1;
  vtxInfo1.constraintVertex = vtx1Constr;
  vtxInfo1.oldPosition = vtxInfo1.linPoint;
  vtxInfo1.seedPosition = vtxInfo1.linPoint;

  for (const auto& trk : params1) {
    vtxInfo1.trackLinks.push_back(&trk);
    state.tracksAtVerticesMap.insert(
        std::make_pair(std::make_pair(&trk, &vtx1),
                       TrackAtVertex<BoundTrackParameters>(1.5, trk, &trk)));
  }

  // Prepare second vertex
  Vector3 vtxPos2(0.3_mm, -0.2_mm, -4.8_mm);
  Vertex<BoundTrackParameters> vtx2(vtxPos2);

  // Add to vertex list
  vtxList.push_back(&vtx2);

  // The constraint vtx for vtx2
  Vertex<BoundTrackParameters> vtx2Constr(vtxPos2);
  vtx2Constr.setFullCovariance(covConstr);
  vtx2Constr.setFitQuality(0, -3);

  // Prepare vtx info for fitter
  VertexInfo<BoundTrackParameters> vtxInfo2;
  vtxInfo2.linPoint.setZero();
  vtxInfo2.linPoint.head<3>() = vtxPos2;
  vtxInfo2.constraintVertex = vtx2Constr;
  vtxInfo2.oldPosition = vtxInfo2.linPoint;
  vtxInfo2.seedPosition = vtxInfo2.linPoint;

  for (const auto& trk : params2) {
    vtxInfo2.trackLinks.push_back(&trk);
    state.tracksAtVerticesMap.insert(
        std::make_pair(std::make_pair(&trk, &vtx2),
                       TrackAtVertex<BoundTrackParameters>(1.5, trk, &trk)));
  }

  state.vtxInfoMap[&vtx1] = std::move(vtxInfo1);
  state.vtxInfoMap[&vtx2] = std::move(vtxInfo2);

  state.addVertexToMultiMap(vtx1);
  state.addVertexToMultiMap(vtx2);

  // Fit vertices
  fitter.fit(state, vtxList, linearizer, vertexingOptions);

  auto vtx1Pos = state.vertexCollection.at(0)->position();
  auto vtx1Cov = state.vertexCollection.at(0)->covariance();
  // auto vtx1Trks = state.vertexCollection.at(0)->tracks();
  auto vtx1FQ = state.vertexCollection.at(0)->fitQuality();

  auto vtx2Pos = state.vertexCollection.at(1)->position();
  auto vtx2Cov = state.vertexCollection.at(1)->covariance();
  // auto vtx2Trks = state.vertexCollection.at(1)->tracks();
  auto vtx2FQ = state.vertexCollection.at(1)->fitQuality();

  if (debugMode) {
    // Vertex 1
    std::cout << "Vertex 1, position: " << vtx1Pos << std::endl;
    std::cout << "Vertex 1, covariance: " << vtx1Cov << std::endl;
    // for (auto t : vtx1Trks) {
    //   std::cout << "\tTrackWeight:" << t.trackWeight << std::endl;
    // }
    std::cout << "Vertex 1, chi2: " << vtx1FQ.first << std::endl;
    std::cout << "Vertex 1, ndf: " << vtx1FQ.second << std::endl;

    // Vertex 2
    std::cout << "Vertex 2, position: " << vtx2Pos << std::endl;
    std::cout << "Vertex 2, covariance: " << vtx2Cov << std::endl;
    // for (auto t : vtx2Trks) {
    //   std::cout << "\tTrackWeight:" << t.trackWeight << std::endl;
    // }
    std::cout << "Vertex 2, chi2: " << vtx2FQ.first << std::endl;
    std::cout << "Vertex 2, ndf: " << vtx2FQ.second << std::endl;
  }

  // Expected values from Athena implementation
  // Vertex 1
  const Vector3 expVtx1Pos(0.077_mm, -0.189_mm, 2.924_mm);

  // Helper matrix to create const expVtx1Cov below
  SymMatrix3 expVtx1Cov;
  expVtx1Cov << 0.329, 0.016, -0.035, 0.016, 0.250, 0.085, -0.035, 0.085, 0.242;

  ActsVector<6> expVtx1TrkWeights;
  expVtx1TrkWeights << 0.8128, 0.7994, 0.8164, 0.8165, 0.8165, 0.8119;
  const double expVtx1chi2 = 0.9812;
  const double expVtx1ndf = 6.7474;

  // Vertex 2
  const Vector3 expVtx2Pos(-0.443_mm, -0.044_mm, -4.829_mm);
  // Helper matrix to create const expVtx2Cov below
  SymMatrix3 expVtx2Cov;
  expVtx2Cov << 1.088, 0.028, -0.066, 0.028, 0.643, 0.073, -0.066, 0.073, 0.435;

  const Vector3 expVtx2TrkWeights(0.8172, 0.8150, 0.8137);
  const double expVtx2chi2 = 0.2114;
  const double expVtx2ndf = 1.8920;

  // Compare the results
  // Vertex 1
  CHECK_CLOSE_ABS(vtx1Pos, expVtx1Pos, 0.001_mm);
  CHECK_CLOSE_ABS(vtx1Cov, expVtx1Cov, 0.001_mm);
  for (int i = 0; i < expVtx1TrkWeights.size(); i++) {
    // CHECK_CLOSE_ABS(vtx1Trks[i].trackWeight, expVtx1TrkWeights[i], 0.001);
  }
  CHECK_CLOSE_ABS(vtx1FQ.first, expVtx1chi2, 0.001);
  CHECK_CLOSE_ABS(vtx1FQ.second, expVtx1ndf, 0.001);

  // Vertex 2
  CHECK_CLOSE_ABS(vtx2Pos, expVtx2Pos, 0.001_mm);
  CHECK_CLOSE_ABS(vtx2Cov, expVtx2Cov, 0.001_mm);
  for (int i = 0; i < expVtx2TrkWeights.size(); i++) {
    // CHECK_CLOSE_ABS(vtx2Trks[i].trackWeight, expVtx2TrkWeights[i], 0.001);
  }
  CHECK_CLOSE_ABS(vtx2FQ.first, expVtx2chi2, 0.001);
  CHECK_CLOSE_ABS(vtx2FQ.second, expVtx2ndf, 0.001);
}

}  // namespace Test
}  // namespace Acts
