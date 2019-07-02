// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE MultiAdaptiveVertexFitter Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Vertexing/MultiAdaptiveVertexFitter.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

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
std::uniform_real_distribution<> pTDist(1. * units::_GeV, 30. * units::_GeV);
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

/// @brief Unit test for MultiAdaptiveVertexFitter
///
BOOST_AUTO_TEST_CASE(multi_adaptive_vertex_fitter_test) {
  bool debugMode = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

  // Set up Eigenstepper
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  Propagator<EigenStepper<ConstantBField>> propagator(stepper);

  VertexFitterOptions<BoundParameters> fitterOptions(tgContext, mfContext);

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>::Config
      config(bField, propagator);

  // Test smoothing
  config.doSmoothing = true;

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>
      fitter(config);

  MultiAdaptiveVertexFitter<ConstantBField, BoundParameters,
                            Propagator<EigenStepper<ConstantBField>>>::State
      state;

  // Create positions of three vertices, two of which (1 and 2) are
  // close to one another and will share a common track later
  Vector3D vtxPos1(-0.15 * units::_mm, -0.1 * units::_mm, -1.5 * units::_mm);
  Vector3D vtxPos2(-0.1 * units::_mm, -0.15 * units::_mm, -3. * units::_mm);
  Vector3D vtxPos3(0.2 * units::_mm, 0.2 * units::_mm, 10. * units::_mm);

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
    std::unique_ptr<Covariance> covMat = std::make_unique<Covariance>();

    (*covMat) << resD0 * resD0, 0., 0., 0., 0., 0., 0., resZ0 * resZ0, 0., 0.,
        0., 0., 0., 0., resPh * resPh, 0., 0., 0., 0., 0., 0., resTh * resTh,
        0., 0., 0., 0., 0., 0., resQp * resQp, 0., 0., 0., 0., 0., 0., 1.;

    // Index of current vertex
    int vtxIdx = (int)(iTrack / nTracksPerVtx);

    // Construct random track parameters
    TrackParametersBase::ParVector_t paramVec;
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

  auto res1 = fitter.addVertexToFit(state, vtxList[0], fitterOptions);

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
                  1 * units::_mm);
  CHECK_CLOSE_ABS(vtxList[1].fullPosition(), seedListCopy[1].fullPosition(),
                  1 * units::_mm);

  auto res2 = fitter.addVertexToFit(state, vtxList[2], fitterOptions);

  BOOST_CHECK(res2.ok());

  // Now also the third vertex should have been modified and fitted
  BOOST_CHECK_NE(vtxList[2].fullPosition(), seedListCopy[2].fullPosition());
  CHECK_CLOSE_ABS(vtxList[2].fullPosition(), seedListCopy[2].fullPosition(),
                  1 * units::_mm);

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

}  // namespace Test
}  // namespace Acts
