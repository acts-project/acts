// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE IterativeVertexFinder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/TrackToVertexIPEstimator.hpp"

namespace bdata = boost::unit_test::data;

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
std::uniform_real_distribution<> resQoPDist(-0.01, 0.01);
// Number of vertices per test event distribution
std::uniform_int_distribution<> nVertexDist(1, 6);
// Number of tracks per vertex distribution
std::uniform_int_distribution<> nTracksDist(5, 15);

///
/// @brief Unit test for IterativeVertexFinder
///
BOOST_AUTO_TEST_CASE(iterative_finder_test) {
  bool debug = false;

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);

  // Number of test events
  unsigned int nEvents = 5;  // = nTest

  for (unsigned int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // Set up constant B-Field
    ConstantBField bField(Vector3D(0., 0., 1.) * units::_T);

    // Set up Eigenstepper
    EigenStepper<ConstantBField> stepper(bField);

    // Set up propagator with void navigator
    Propagator<EigenStepper<ConstantBField>> propagator(stepper);

    typedef FullBilloirVertexFitter<ConstantBField, BoundParameters,
                                    Propagator<EigenStepper<ConstantBField>>>
        BilloirFitter;

    // Set up Billoir Vertex Fitter
    BilloirFitter::Config vertexFitterCfg(bField, propagator);

    std::unique_ptr<BilloirFitter> bFitterPtr

        = std::make_unique<BilloirFitter>(vertexFitterCfg);

    // Vertex Finder
    typedef IterativeVertexFinder<ConstantBField, BoundParameters,
                                  Propagator<EigenStepper<ConstantBField>>>
        VertexFinder;
    VertexFinder::Config cfg(bField, std::move(bFitterPtr), propagator);
    cfg.reassignTracksAfterFirstFit = true;

    VertexFinder finder(cfg);

    // Vector to be filled with all tracks in current event
    std::vector<BoundParameters> tracks;

    // Vector to be filled with truth vertices for later comparison
    std::vector<Vertex<BoundParameters>> trueVertices;

    // start creating event with nVertices vertices
    unsigned int nVertices = nVertexDist(gen);
    for (unsigned int iVertex = 0; iVertex < nVertices; ++iVertex) {
      // Number of tracks
      unsigned int nTracks = nTracksDist(gen);

      if (debug) {
        std::cout << "Event " << iEvent << ", Vertex " << iVertex << "/"
                  << nVertices << " with " << nTracks << " tracks."
                  << std::endl;
      }
      // Create perigee surface
      std::shared_ptr<PerigeeSurface> perigeeSurface =
          Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      // True vertex
      Vertex<BoundParameters> trueV(Vector3D(x, y, z));
      std::vector<TrackAtVertex<BoundParameters>> tracksAtTrueVtx;

      // Calculate d0 and z0 corresponding to vertex position
      double d0_v = sqrt(x * x + y * y);
      double z0_v = z;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        TrackParametersBase::ParVector_t paramVec;
        double z0track = z0_v + z0Dist(gen);
        paramVec << d0_v + d0Dist(gen), z0track, phiDist(gen), thetaDist(gen),
            q / pTDist(gen), 0.;

        // Fill vector of track objects with simple covariance matrix
        std::unique_ptr<Covariance> covMat = std::make_unique<Covariance>();

        // Resolutions
        double res_d0 = resIPDist(gen);
        double res_z0 = resIPDist(gen);
        double res_ph = resAngDist(gen);
        double res_th = resAngDist(gen);
        double res_qp = resQoPDist(gen);

        (*covMat) << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0,
            0., 0., 0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
            res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0.,
            0., 0., 0., 0., 1.;
        auto params = BoundParameters(tgContext, std::move(covMat), paramVec,
                                      perigeeSurface);
        tracks.push_back(params);

        TrackAtVertex<BoundParameters> trAtVt(0., params, params);
        tracksAtTrueVtx.push_back(trAtVt);
      }

      trueV.setTracksAtVertex(tracksAtTrueVtx);
      trueVertices.push_back(trueV);

    }  // end loop over vertices

    // shuffle list of tracks
    std::shuffle(std::begin(tracks), std::end(tracks), gen);

    VertexFinderOptions<BoundParameters> vFinderOptions(tgContext, mfContext);

    // find vertices
    auto res = finder.find(tracks, vFinderOptions);

    BOOST_CHECK(res.ok());

    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    // check if same amount of vertices has been found with tolerance of 2
    CHECK_CLOSE_ABS(vertexCollection.size(), nVertices, 2);

    if (debug) {
      std::cout << "########## RESULT: ########## Event " << iEvent
                << std::endl;
      std::cout << "Number of true vertices: " << nVertices << std::endl;
      std::cout << "Number of reco vertices: " << vertexCollection.size()
                << std::endl;

      int count = 1;
      std::cout << "----- True vertices -----" << std::endl;
      for (const auto& vertex : trueVertices) {
        Vector3D pos = vertex.position();
        std::cout << count << ". True Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
      std::cout << "----- Reco vertices -----" << std::endl;
      count = 1;
      for (const auto& vertex : vertexCollection) {
        Vector3D pos = vertex.position();
        std::cout << count << ". Reco Vertex:\t Position:"
                  << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
                  << std::endl;
        std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
                  << std::endl;
        count++;
      }
    }

    // Check if all vertices have been found with close z-values
    bool allVerticesFound = true;
    for (const auto& trueVertex : trueVertices) {
      SpacePointVector truePos = trueVertex.fullPosition();
      bool currentVertexFound = false;
      for (const auto& recoVertex : vertexCollection) {
        SpacePointVector recoPos = recoVertex.fullPosition();
        // check only for close z distance
        double zDistance = std::abs(truePos[eZ] - recoPos[eZ]);
        if (zDistance < 2 * units::_mm) {
          currentVertexFound = true;
        }
      }
      if (!currentVertexFound) {
        allVerticesFound = false;
      }
    }

    // check if found vertices have compatible z values
    BOOST_TEST(allVerticesFound);
  }
}

}  // namespace Test
}  // namespace Acts
