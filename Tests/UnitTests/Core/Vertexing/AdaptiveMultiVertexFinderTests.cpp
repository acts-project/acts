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

#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include "AMVFTestData.ipp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_CASE(adaptive_multi_vertex_finder_test) {
  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);
  PropagatorOptions<> pOptions(tgContext, mfContext);
  pOptions.direction = backward;

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  // IP 3D Estimator
  using IPEstimator = ImpactPoint3dEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ip3dEstCfg(bField, propagator, pOptions, false);
  IPEstimator ip3dEst(ip3dEstCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundParameters, Linearizer>;

  Fitter::Config fitterCfg(ip3dEst);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator, pOptions);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = false;

  Fitter fitter(fitterCfg);

  using SeedFinder = TrackDensityVertexFinder<Fitter, GaussianTrackDensity>;

  SeedFinder seedFinder;

  using IPEstimater = TrackToVertexIPEstimator<BoundParameters, Propagator>;

  IPEstimater::Config ipEstCfg(propagator, pOptions);

  // Create TrackToVertexIPEstimator
  IPEstimater ipEst(ipEstCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), std::move(seedFinder),
                              std::move(ipEst), std::move(linearizer));

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  Finder finder(finderConfig);

  auto tracks = getAthenaTracks();

  if (debugMode) {
    std::cout << "Number of tracks in event: " << tracks.size() << std::endl;
    int maxCout = 10;
    int count = 0;
    for (const auto& trk : tracks) {
      std::cout << count << ". track: " << std::endl;
      std::cout << "params: " << trk << std::endl;
      count++;
      if (count == maxCout) {
        break;
      }
    }
  }

  std::vector<const BoundParameters*> tracksPtr;
  for (const auto& trk : tracks) {
    tracksPtr.push_back(&trk);
  }

  VertexFinderOptions<BoundParameters> finderOptions(tgContext, mfContext);

  Vector3D constraintPos{0._mm, 0._mm, 0_mm};
  ActsSymMatrixD<3> constraintCov;
  constraintCov << 0.000196000008145347238, 0, 0, 0, 0.000196000008145347238, 0,
      0, 0, 2809;

  Vertex<BoundParameters> constraintVtx;
  constraintVtx.setPosition(constraintPos);
  constraintVtx.setCovariance(constraintCov);

  finderOptions.vertexConstraint = constraintVtx;

  auto findResult = finder.find(tracksPtr, finderOptions);

  if (!findResult.ok()) {
    std::cout << findResult.error().message() << std::endl;
  }

  BOOST_CHECK(findResult.ok());

  std::vector<Vertex<BoundParameters>> allVertices = *findResult;

  if (debugMode) {
    std::cout << "Number of vertices reconstructed: " << allVertices.size()
              << std::endl;

    int count = 0;
    for (const auto& vtx : allVertices) {
      count++;
      std::cout << count << ". Vertex at position: " << vtx.position()[0]
                << ", " << vtx.position()[1] << ", " << vtx.position()[2]
                << std::endl;
      std::cout << count << ". Vertex with cov: " << vtx.covariance()
                << std::endl;
      std::cout << "\t with n tracks: " << vtx.tracks().size() << std::endl;
    }
  }

  // Test expected outcomes from athena implementation
  // Number of reconstructed vertices
  const int expNRecoVertices = 15;

  // First vertex
  const Vector3D expVtx1Pos(-0.0067_mm, 0.0060_mm, -6.0709_mm);
  ActsSymMatrixD<3> expVtx1Cov;
  expVtx1Cov << 0.000, 1.e-05, -8.e-05, 1.e-05, 0.000, -8.e-05, -8.e-05,
      -8.e-05, 0.002;
  std::vector<double> expVtx1TrkWeights{0.9796, 0.0334, 0.9884, 0.9697};
  std::vector<double> expVtx1TrkComp{1.2542, 15.7317, 0.1144, 2.067};
  std::vector<double> expVtx1TrkChi2{0, 0, 0, 0};

  // Last vertex
  const Vector3D expVtx15Pos(0.00264_mm, -0.0072_mm, -39.8197_mm);
  ActsSymMatrixD<3> expVtx15Cov;
  expVtx15Cov << 0.000, 1.e-06, 0.000, 1.e-06, 0.000, -6.e-05, 0.000, -6.e-05,
      0.014;
  std::vector<double> expVtx15TrkWeights{0.0048, 0.0005, 0.0236, 0.8481,
                                         0.8924};
  std::vector<double> expVtx15TrkComp{19.6561, 24.1389, 16.4425, 5.5604,
                                      4.7683};
  std::vector<double> expVtx15TrkChi2{0, 0, 0, 0};

  // Vertex z positions of all found vertices
  const std::vector<double> expAllVtxZPos{
      -6.070_mm,   -12.0605_mm, -15.1093_mm, -27.6569_mm, -22.1054_mm,
      -45.7010_mm, -5.0622_mm,  -26.5496_mm, -28.9597_mm, -37.7430_mm,
      5.4828_mm,   -47.8939_mm, 2.5777_mm,   -0.2656_mm,  -39.8197_mm};

  // Number of tracks of all vertices
  const std::vector<int> expAllNTracks{4, 2, 3, 14, 5, 9, 8, 17,
                                       7, 2, 2, 4,  2, 7, 5};

  BOOST_CHECK_EQUAL(allVertices.size(), expNRecoVertices);

  int count = 0;
  for (const auto& vtx : allVertices) {
    // Check vertex z positions
    CHECK_CLOSE_ABS(vtx.position()[2], expAllVtxZPos[count], 0.003_mm);
    // Check number of tracks
    BOOST_CHECK_EQUAL(vtx.tracks().size(), expAllNTracks[count]);

    // Check vertex 1 thoroughly
    if (count == 0) {
      CHECK_CLOSE_ABS(vtx.position(), expVtx1Pos, 0.001_mm);
      CHECK_CLOSE_ABS(vtx.covariance(), expVtx1Cov, 0.001_mm);
      int trkCount = 0;
      for (const auto& trk : vtx.tracks()) {
        CHECK_CLOSE_ABS(trk.trackWeight, expVtx1TrkWeights[trkCount], 0.01);
        CHECK_CLOSE_ABS(trk.vertexCompatibility, expVtx1TrkComp[trkCount],
                        0.15);
        // CHECK_CLOSE_ABS(trk.chi2Track, expVtx1TrkChi2[trkCount], 0.001);
        trkCount++;
      }
    }

    // Check vertex 15 thoroughly
    if (count == 14) {
      CHECK_CLOSE_ABS(vtx.position(), expVtx15Pos, 0.001_mm);
      CHECK_CLOSE_ABS(vtx.covariance(), expVtx15Cov, 0.001_mm);
      int trkCount = 0;
      for (const auto& trk : vtx.tracks()) {
        CHECK_CLOSE_ABS(trk.trackWeight, expVtx15TrkWeights[trkCount], 0.01);
        CHECK_CLOSE_ABS(trk.vertexCompatibility, expVtx15TrkComp[trkCount],
                        0.15);
        // CHECK_CLOSE_ABS(trk.chi2Track, expVtx15TrkChi2[trkCount], 0.001);
        trkCount++;
      }
    }

    count++;
  }
}

}  // namespace Test
}  // namespace Acts
